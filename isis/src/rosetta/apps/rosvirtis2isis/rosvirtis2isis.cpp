#include "Isis.h"

#include <iostream>
#include <istream>
#include <sstream>
#include <cmath>

#include <QString>
#include <QFileInfo>
#include <QByteArray>
#include <QTextStream>
#include <QStringList>

#include "FileName.h"
#include "ImportPdsTable.h"
#include "ProcessImportPds.h"
#include "Table.h"
#include "UserInterface.h"
#include "PolynomialUnivariate.h"
#include "VirtisHK.h"

using namespace std;
using namespace Isis;

// helper functions
QByteArray pvlFix(QString fileName); 
int word(int byte1, int byte2);
double translateScet(int word1, int word2, int word3);
bool isValid(int word);

void IsisMain ()
{
  ProcessImportPds p;
  UserInterface &ui = Application::GetUserInterface();

  FileName inFile = ui.GetFileName("FROM");

  Pvl pdsLabel; 
  try {
    pdsLabel.read(inFile.expanded());
  } 
  catch (IException &e) {
    // Try to fix the PVL before reading it in
    QByteArray pvlData = pvlFix(inFile.expanded()); 
    QTextStream pvlTextStream(&pvlData);
    istringstream pvlStream(pvlTextStream.readAll().toStdString());

    try {
      pvlStream >> pdsLabel;
    }
    catch(IException &e) {
      QString msg = "Input file [" + inFile.expanded() +
                 "] is not a valid PVL file.";
      throw IException(e, IException::User, msg, _FILEINFO_);
    }
  }
  
  p.SetPdsFile(pdsLabel, inFile.expanded()); 
  p.SetOrganization(Isis::ProcessImport::BIP);

  // NULL pixels are set to 65535 in the input QUB
  p.SetNull(65535, 65535); 

  Cube *outcube = p.SetOutputCube ("TO");

  // Is this a correctly-formatted Rosetta VIRTIS-M file? (VIRTIS-H is not currently supported)
  QString instid, missid, channelid; 

  try {
    // ROSETTA:CHANNEL_ID will be IR or VIS
    instid = (QString) pdsLabel.findKeyword("INSTRUMENT_ID"); 
    missid = (QString) pdsLabel.findKeyword("MISSION_ID");
    channelid = (QString) pdsLabel.findKeyword("ROSETTA:CHANNEL_ID");
  }
  catch (IException &e) {
    QString msg = "Unable to read [INSTRUMENT_ID] or [MISSION_ID] from input file [" +
                 inFile.expanded() + "]";
    throw IException(e, IException::Io,msg, _FILEINFO_);
  }

  instid = instid.simplified().trimmed();
  missid = missid.simplified().trimmed();
  channelid = channelid.simplified().trimmed(); 

  if (missid != "ROSETTA" && instid != "VIRTIS") {
    QString msg = "Input file [" + inFile.expanded() + "] does not appear to be a " +
                 "DAWN Visual and InfraRed Mapping Spectrometer (VIR) EDR or RDR file.";
    throw IException(IException::Unknown, msg, _FILEINFO_);
  }

  // Override default DataTrailerBytes constructed from PDS header
  // Will this number ever change? Where did this # come from? 
  p.SetDataTrailerBytes(864); 

  p.StartProcess();

  // Retrieve HK settings file and read in HK values. 
  QList<VirtisHK> hk; 

  QFile hkFile("./assets/virtis_housekeeping.txt");

  if(!hkFile.open(QIODevice::ReadOnly)) {
    // Appropriate error
  }

  QTextStream in(&hkFile);

  while(!in.atEnd()) {
      QString line = in.readLine();    
      QStringList fields = line.split(",");    
      hk.append(VirtisHK(fields[0], fields[1], fields[2], fields[3], fields[4]));
  }

  hkFile.close();

  // Construct HK (housekeeping) table
  TableRecord rec; 

  QList<TableField> tableFields; 

  for (int i=0; i < hk.size(); i++) {
    tableFields.append(hk[i].tableField()); 
  }

  for (int i=0; i < tableFields.size(); i++) {
    rec += tableFields[i];
  }

  Table table("VIRTISHouseKeeping", rec);
                                                  
  // VIRTIS-M (VIS and IR) Equations
  // These are adapted from the VIRTIS IDL processing pipeline
  // and pg. 66-67 of the VIRTIS-EAICD
  QList<std::vector<double> > equationList; 
  for (int i=0; i < hk.size(); i++) {
    equationList.append(hk[i].coefficients()); 
  }

  QList<PolynomialUnivariate> equations;

  for (int s=0; s < equationList.size(); s++) {
    equations.append(PolynomialUnivariate(2, equationList[s])); 
  }
  
  // Populate the Housekeeping table
  // 
  // There are 3 categories of VIRTIS HK Values, in terms of converting from input byte to output 
  // value:
  // 
  // (1) SCET (many-to-one)
  // (2) Physical Quantities (one-to-one)
  // (3) Flags or Statistics (one-to-many)
  // 
  // SCET values are made up of 3 VIRTIS HK 2-byte words. The output value can be calculated
  // using the translateScet helper function.
  // 
  // Physical values are made up of 1 VIRTIS HK 2-byte word, which is converted to a physical value
  // using an equation specified as a series of coefficients in the associated "assets" file. 
  //
  // For Flags or Statistics Values, 1 VIRTIS HK 2-byte word is associated with several (a varaible
  // number of) output Flag or Statistics values. These are all treated as special cases.
  //
  // Additionally, Sine and Cosine HK Values need to be pre-processed before conversion, but are 
  // otherwise treated as a normal "Physical Quantity" HK. 
  //
  std::vector< char * > hkData = p.DataTrailer(); 
  for (unsigned int i=0; i < hkData.size() ; i++) {
    const char *hk = hkData.at(i);

    // Each data trailer can contain multiple 82-word records, but we only need 1 / line
    int start = 0;
    int tableNum = 0;     

    // Loop through each 82-word record
    // Each k is a byteNumber
    for (int k=0; k<82*2; k=k+2) {
      int temp = 0;

      // Convert non-SCET records (1 two-byte word each) using the appropriate equations
      if (k !=0 && k!=14 && k!=38 && k!=58 && k!=116) {
        temp = word(hk[start + k], hk[start+k+1]);
        if (temp != 65535) {

          // If Sine or Cosine, pre-process before sending to conversion.
          if (tableNum == 63) { // SIN
            int HK_bit = (int) (temp) & 4095;
            int HK_sign = (int) (temp/4096) & 1; 
            temp = (int) HK_bit * HK_sign;
          }
          if (tableNum == 64) { // COS
            temp = (int) (temp) & 4095;
          }

          // Specical one-to-many cases (Flags or Statistics) 
          if (tableNum == 2) { // # of Subslices / first seial num 2-3
            rec[tableNum] = hk[start+k]*1.0; 
            rec[tableNum+1] = hk[start+k+1]*1.0;
            tableNum++;
          } else if (tableNum == 4) { // Data Type 4-9
            rec[tableNum] = (int)(temp/-32768) & 1; 
            rec[tableNum+1] = (int)(temp/16384) & 1;
            rec[tableNum+2] = (int)(temp/8192) & 1;
            rec[tableNum+3] = (int)(temp/1024) & 7;
            rec[tableNum+4] = (int)(temp/256) & 3;
            rec[tableNum+5] = (int)(temp/1) & 255;
            tableNum+=5;
          } else if (tableNum == 12) { // V_MODE 12-14
            rec[tableNum] = 1.0* ((int)(temp/4096) & 15);
            rec[tableNum+1] = 1.0* ((int)(temp/64) & 63);
            rec[tableNum+2] = 1.0* ((int)(temp/1) & 63);
            tableNum+=2;
          } else if (tableNum == 15) { //ME_PWR_STAT 15 - 21
            rec[tableNum] = 1.0* ((int)(temp/1) & 1);
            rec[tableNum+1] = 1.0* ((int)(temp/2) & 1);
            rec[tableNum+2] = 1.0* ((int)(temp/4) & 1);
            rec[tableNum+3] = 1.0* ((int)(temp/8) & 1);
            rec[tableNum+4] = 1.0* ((int)(temp/16) & 1);
            rec[tableNum+5] = 1.0* ((int)(temp/32) & 1);
            rec[tableNum+6] = 1.0* ((int)(temp/-32786) & 1);
            tableNum+=6;
          } else if (tableNum == 30){  // M_ECA_STAT 30-31
            rec[tableNum] = 1.0* ((int)(temp/1) & 1);
            rec[tableNum+1] = 1.0* ((int)(temp/256) & 1);
            tableNum++;
          } else if (tableNum == 32) { // M_COOL_STAT 32-34
            rec[tableNum] = 1.0* ((int)(temp/1) & 1);
            rec[tableNum+1] = 1.0* ((int)(temp/16) & 1);
            rec[tableNum+2] = 1.0* ((int)(temp/256) & 1);
            tableNum+=2;
          } 

          else if (tableNum == 65) { // M_VIS_FLAG
            rec[tableNum] = 1.0* ((int)(temp/1) & 1);
            rec[tableNum+1] = 1.0* ((int)(temp/2) & 1);
            rec[tableNum+2] = 1.0* ((int)(temp/4) & 1);
            rec[tableNum+3] = 1.0* ((int)(temp/8) & 1);
            rec[tableNum+4] = 1.0* ((int)(temp/16) & 1);
            rec[tableNum+5] = 1.0* ((int)(temp/256) & 1);
            tableNum+=5;
          } 
          else if (tableNum == 91) { //M_IR_LAMP_SHUTTER
            double lamp1 = 1.0* ((int)(temp/1) & 15);
            rec[tableNum] = equations[tableNum].Evaluate(lamp1);
            rec[tableNum+1] = 1.0* ((int)(temp/16) & 1);
            double lamp2 = 1.0* ((int)(temp/256) & 15);
            rec[tableNum+2] = equations[tableNum+1].Evaluate(lamp2);
            rec[tableNum+3] = 1.0* ((int)(temp/4096) & 1);
            tableNum+=3;
          } 
          else if (tableNum == 95) { // M_IR_FLAG
            rec[tableNum] = 1.0* ((int)(temp/1) & 1);
            rec[tableNum+1] = 1.0* ((int)(temp/2) & 1);
            rec[tableNum+2] = 1.0* ((int)(temp/4) & 1);
            rec[tableNum+3] = 1.0* ((int)(temp/8) & 1);
            rec[tableNum+4] = 1.0* ((int)(temp/16) & 1);
            rec[tableNum+5] = 1.0* ((int)(temp/32) & 1);
            rec[tableNum+6] = 1.0* ((int)(temp/64) & 1);
            rec[tableNum+7] = 1.0* ((int)(temp/512) & 1);
            rec[tableNum+8] = 1.0* ((int)(temp/4096) & 1);
            rec[tableNum+9] = 1.0* ((int)(temp/8192) & 1);
            rec[tableNum+10] = 1.0* ((int)(temp/16384) & 1);
            tableNum+=10;
          }
          else {
             // Convert a physical quantity to its output value (1 word -> 1 output physical value)
             rec[tableNum] = equations[tableNum].Evaluate(temp * 1.0); 
          }
        }
        else {
          rec[tableNum] = 65535.0; // HK Data is Invalid
        }
      }
      else { 
        // Convert SCET records (3 words -> one output SCET) 
        int word1 = word(hk[start+k], hk[start+k+1]);
        int word2 = word(hk[start+k+2], hk[start+k+3]);
        int word3 = word(hk[start+k+4], hk[start+k+5]);

        int result;

        // If any of the words comprising the SCET are invalid, the whole thing is invalid. 
        if (isValid(word1) && isValid(word2) && isValid(word3)) {
          result = translateScet(word1, word2, word3); 
        }
        else {
          result = 65535; 
        }
        
        // If we don't have a valid SCET, the whole line of HK data is not valid, so we skip it.
        if (result == 0 || result == 65535) {
          break; 
        }
        else{
          rec[tableNum] = result*1.0; 
        }
        // We used 3 words
        k=k+4; 
      }
      tableNum++;
    }
    table += rec; 
  }

  outcube->write(table);

  // Get the directory where the Rosetta translation tables are.
  PvlGroup dataDir (Preference::Preferences().findGroup("DataDirectory"));
  QString transDir = (QString) dataDir["Rosetta"] + "/translations/";

  // Create a PVL to store the translated labels in
  Pvl outLabel;

  // Translate the Archive group
  FileName transFile = transDir + "virtisArchive.trn";
  PvlTranslationManager archiveXlater (pdsLabel, transFile.expanded());
  archiveXlater.Auto(outLabel);

  // Translate the Instrument group
  transFile = transDir + "virtisInstruments.trn";
  PvlTranslationManager instrumentXlater (pdsLabel, transFile.expanded());
  instrumentXlater.Auto(outLabel);

  // Write the Archive and Instrument groups to the output cube label
  outcube->putGroup(outLabel.findGroup("Archive", Pvl::Traverse));
  outcube->putGroup(outLabel.findGroup("Instrument", Pvl::Traverse));

  // Add correct kernels and output Kernel group for VIRTIS_M IR vs. VIS. 
  // VIRTIS_H is not supported.
  PvlGroup kerns("Kernels");
  if (channelid == "VIRTIS_M_IR") {
    kerns += PvlKeyword("NaifFrameCode", toString(-226212));
  }
  else if (channelid == "VIRTIS_M_VIS") {
    kerns += PvlKeyword("NaifFrameCode", toString(-226214));
  } 
  else {
    QString msg = "Input file [" + inFile.expanded() + "] has an invalid " +
                 "InstrumentId.";
    throw IException(IException::Unknown, msg, _FILEINFO_);
  }
  outcube->putGroup(kerns);

  p.EndProcess ();
}


/**
 * Convert 2 bytes into a 2-byte word
 * 
 * @param byte1 first byte (MSB)
 * @param byte2 second byte
 * 
 * @return int 2-word byte
 */
int word(int byte1, int byte2)
{
  int result = 65535;
  if (!(((int) byte1 == -1) && ((int) byte2 == -1))) {
    result = 0; 
    result += ((int)byte1) << 8; 
    result += ((int) byte2);
  }
  return result;  
}


/**
 * Translate the three constituent VIRTIS HK words into a 3-word SCET 
 * (SpaceCraft Event Time) value.  
 * 
 * @param word1 first word (most-significant word)
 * @param word2 second word
 * @param word3 third word (least-significant word)
 * 
 * @return int calculated SCET value
 */
double translateScet(int word1, int word2, int word3)
{
  return (unsigned long) (word1* pow(2.0,16.0)) + (unsigned long) word2 + (double) (word3)/pow(2.0,16.0);
}



/**
 * True if valid, false otherwise
 * 
 * @author kberry (5/25/2017)
 * 
 * @param word The VIRTIS housekeeping word
 * 
 * @return bool True if valid, else False
 */
bool isValid(int word){
  return word != 65535; 
}


/**
 * @brief   Fixes incorrectly formatted PDS headers for VIRTIS data in draft data area 
 *          as of 04/2016. Unnecessary for non-draft data.
 *  
 * @param   fileName The full path to the VIRTIS input cube
 * @return  A QByteArray containing a valid, fixed PDS header. This is fed to
 *          a QTextStream which is then fed into a Pvl object.
 *
 */
QByteArray pvlFix(QString fileName){
  QByteArray null;
  QFile pvlFile;

  pvlFile.setFileName(fileName);

  if (!pvlFile.open(QFile::ReadOnly|QIODevice::Text)) {
    return null;
    //change to except
  }

  //Read the Pvl file into a byte array
  QByteArray fileData = pvlFile.readAll();
  QByteArray pvlData;

  QString pvlEnd("QUBE\n\nEND");
  int ix = fileData.lastIndexOf(pvlEnd);

  pvlData = fileData.left(ix+pvlEnd.size());

  QList<QByteArray> pvlLines = pvlData.split('\n');

  bool withinSpice = false; 

  QByteArray newPvlData; 

  for (int i=0; i< pvlLines.size(); i++) {
    QString local(pvlLines[i]);
    if (local.contains("SOFTWARE_VERSION_ID")) {
      pvlLines[i].append(","); 
    } 
    else if (local.contains("SPICE_FILE_NAME") || withinSpice) {
      if (local.contains(")")) {
        withinSpice = false; 
      }
      else {
        pvlLines[i].append(","); 
        withinSpice = true;
      }
    }
    newPvlData.append(pvlLines[i].append("\n")); 
  }
    return newPvlData;
}

