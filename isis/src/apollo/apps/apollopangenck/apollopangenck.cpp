/**
 *   Unless noted otherwise, the portions of Isis written by the USGS are
 *   public domain. See individual third-party library and package descriptions
 *   for intellectual property information, user agreements, and related
 *   information.
 *
 *   Although Isis has been used by the USGS, no warranty, expressed or
 *   implied, is made by the USGS as to the accuracy and functioning of such
 *   software and related material nor shall the fact of distribution
 *   constitute any such warranty, and no responsibility is assumed by the
 *   USGS in connection therewith.
 *
 *   For additional information, launch
 *   $ISISROOT/doc//documents/Disclaimers/Disclaimers.html
 *   in a browser or see the Privacy &amp; Disclaimers page on the Isis website,
 *   http://isis.astrogeology.usgs.gov, and the USGS privacy and disclaimers on
 *   http://www.usgs.gov/privacy.html.
 */

// Isis
#include "Isis.h"

// C++ standard libraries if needed
#include <math.h>

//Isis Headers if needed
#include "Angle.h"
#include "ApolloPanTile.h"
#include "CSVReader.h"
#include "FileName.h"
#include "iTime.h"
#include "LinearAlgebra.h"
#include "UserInterface.h"

using namespace Isis;

void IsisMain() {
  UserInterface &ui = Application::GetUserInterface();

  double centerFiducialSample = 119695; // sample in fiducial space
  double focalLength = 609.752;         // focal length from Wu, S.C. and Moore, H.J. 1980
  double pixSizeMM = 0.005;             // Experimental Photography of Lunar Images
                                        // Geological Survey Professional Paper 1046-D

  LinearAlgebra::Matrix tmp(3,3);       // temporary matrix for processing
  LinearAlgebra::Matrix phi0(3,3);      // phi rotation matrix (about y) at image center time
  LinearAlgebra::Matrix omega0(3,3);    // omega rotation matrix (about x) at image center time
  LinearAlgebra::Matrix kappa0(3,3);    // kappa rotatation matrix (about y) at image center time
  LinearAlgebra::Matrix R0(3,3);        // full rotation matrix at image center time
  LinearAlgebra::Matrix sweep(3,3);     // sweep angle matrix at any time
  LinearAlgebra::Matrix R(3,3);         // rotation matrix at any time (transpose(sweep) * R0)

  QString str;
  QStringList outputList;               // list of output strings

  FileName tileFileName;
  ApolloPanTile tile;

  // input files
  FileName pvlFileNamePrefix = ui.GetFileName("FROMPVL");
  const QString supportFileName = ui.GetFileName("FROMREVSUPPORT");

  // csv reader
  CSVReader::CSVAxis row;
  CSVReader reader(supportFileName, true, 1, ',', true, true);

  // strip 4 digit imageNumber from filename and remove any leading zeros
  QString imageNumber = pvlFileNamePrefix.name().right(4);
  imageNumber.remove( QRegExp("^[0]*") );

  // find row in csv file for imageNumber
  for (int i = 0; i < reader.rows(); i++) {
    row = reader.getRow(i);

    if( QString(row[0]) == imageNumber)
      break;
  }

  // get phi, kappa, omega rotation angles from support page csv
  // these angles define the orientation at the center time of the pan image
  // phi is rotation about y; omegexpETa about x; kappa about z
  str = row[11];
  Angle phi = Angle(str.toDouble(), Angle::Units::Degrees);
  str = row[12];
  Angle kappa = Angle(str.toDouble(), Angle::Units::Degrees);
  str = row[13];
  Angle omega = Angle(str.toDouble(), Angle::Units::Degrees);

  // cosines and sines of phi, kappa, omega
  double cosp = cos(phi.radians());
  double sinp = sin(phi.radians());
  double cosk = cos(kappa.radians());
  double sink = sin(kappa.radians());
  double coso = cos(omega.radians());
  double sino = sin(omega.radians());

  // clear matrices
  phi0.clear();
  omega0.clear();
  kappa0.clear();

  // phi matrix at image center time
  phi0(0,0) = cosp;
  phi0(0,2) = sinp;
  phi0(1,1) = 1.0;
  phi0(2,0) = -sinp;
  phi0(2,2) = cosp;

  // omega matrix at image center time
  omega0(0,0) = 1.0;
  omega0(1,1) = coso;
  omega0(1,2) = -sino;
  omega0(2,1) = sino;
  omega0(2,2) = coso;

  // kappa matrix at image center time
  kappa0(0,0) = cosk;
  kappa0(0,1) = -sink;
  kappa0(1,0) = sink;
  kappa0(1,1) = cosk;
  kappa0(2,2) = 1.0;

  // form full rotation matrix at image center time
  // product of phi, omega, kappa matrices
  tmp = prod(phi0,omega0);
  R0 = prod(tmp,kappa0);

  // loop over pvls for image
  // TODO: handle images with less than 8 tiles
  for (int i = 0; i < 8; i++) {
    tileFileName = pvlFileNamePrefix.toString() + "_000" + QString::number(i+1) + ".pvl";

    // read tile pvl
    tile.fromPvl(tileFileName.toString());

    // loop over timing marks in tile
    int numTimingMarks = tile.numberOfTimingMarks();
    for ( int j = 0; j < numTimingMarks; j++) {

      double expSample = tile.expSampleTime(j);       // sample at this timing mark
      iTime etime(tile.etime(j));                     // ephemeris time at this timing mark

      // reference sample of timing mark to center fiducial sample
      double relativeSample = centerFiducialSample - expSample;

      // arc length in mm
      double arcLength = relativeSample * pixSizeMM;

      // compute theta (sweep) angle relative to center in radians
      Angle theta = Angle(arcLength/focalLength, Angle::Units::Radians);

      // create sweep matrix
      sweep.clear();
      sweep(0,0) = 1.0;
      sweep(1,1) = cos(theta.radians());
      sweep(1,2) = -sin(theta.radians());
      sweep(2,1) = sin(theta.radians());
      sweep(2,2) = cos(theta.radians());

      // create full rotation matrix
      // product of transpose of sweep matrix and R0
      R = prod(trans(sweep),R0);

      // create output string
      // ephemeris time plus 9 elements of rotation matrix
      str = QString("%1 %2 %3 %4 %5 %6 %7 %8 %9 %10").arg(etime.UTC()).arg(R(0,0)).arg(R(0,1)).
          arg(R(0,2)).arg(R(1,0)).arg(R(1,1)).arg(R(1,2)).arg(R(2,0)).arg(R(2,1)).arg(R(2,2));

      // add string to QStringList for output
      outputList.append(str);
    }
  }

  // sort list in ascending order
  outputList.sort();

  // loop over strings to remove any with duplicate times
  for (int i = 1; i < outputList.size()-1; i++) {
    if (outputList.at(i-1).left(27) == outputList.at(i).left(27)) {
      outputList.removeAt(i-1);
    }
  }

  // check last two lines for duplicates
  int nLines = outputList.size();
  if (outputList.at(nLines-1).left(27) == outputList.at(nLines-2).left(27)) {
    outputList.removeAt(nLines-1);
  }

  // padding necessary for interpolation at time extents
  // pad beginning of list with two lines
  QString padstr = outputList.first();
  iTime startPadTime(padstr.left(27));
  startPadTime -= 0.1;
  padstr.replace(0,27,startPadTime.UTC());
  outputList.prepend(padstr);
  startPadTime -= 0.1;
  padstr.replace(0,27,startPadTime.UTC());
  outputList.prepend(padstr);

  // pad end of list with two lines
  padstr = outputList.last();
  iTime stopPadTime(padstr.left(27));
  stopPadTime += 0.1;
  padstr.replace(0,27,stopPadTime.UTC());
  outputList.append(padstr);
  stopPadTime += 0.1;
  padstr.replace(0,27,stopPadTime.UTC());
  outputList.append(padstr);

  // create output text file name
  FileName outputFileName = pvlFileNamePrefix.name();
  QString ofname = outputFileName.toString() + "_ck.txt";

  // create output stream
  std::ofstream fpOut(ofname.toLatin1().data(), std::ios::out);
  if (!fpOut) {
    return;
  }

  // output
  for (int i = 0; i < outputList.size(); i++) {
     fpOut << outputList.at(i) << "\n";
  }

  // close output stream
  fpOut.close();

  return;
}
