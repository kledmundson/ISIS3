#include "Isis.h"

#include <iomanip>

#include "Distance.h"
#include "ProcessByLine.h"
#include "Projection.h"
#include "SpecialPixel.h"
#include "LineManager.h"
#include "OriginalLabel.h"
#include "History.h"
#include "Table.h"

using namespace std;
using namespace Isis;

void DoWrap(Buffer &in);
void GetStats(Buffer &in, Buffer &out);

Cube *ocube;

int leftPad;
int rightPad;
int topPad;
int bottomPad;
int inl;
Statistics inCubeStats;
Statistics outCubeStats;

void IsisMain() {
  // We will be using a mosaic technique so get the size of the input file
  ProcessByLine p;
  UserInterface &ui = Application::GetUserInterface();
  Cube *icube = p.SetInputCube("FROM");
  int ins = icube->getSampleCount();
  inl = icube->getLineCount();
  int inb = icube->getBandCount();
  outCubeStats.Reset();

  PvlGroup mapgrp = icube->getLabel()->FindGroup("Mapping", Pvl::Traverse);
  bool hasExtents = false;
  bool isGlobal = false;
  double minLat,maxLat,minLon,maxLon;
  if (mapgrp.HasKeyword("MinimumLatitude") && mapgrp.HasKeyword("MaximumLatitude") &&
      mapgrp.HasKeyword("MinimumLongitude") && mapgrp.HasKeyword("MaximumLongitude")) {
    hasExtents = true;
    minLat = mapgrp["MinimumLatitude"];
    maxLat = mapgrp["MaximumLatitude"];
    minLon = mapgrp["MinimumLongitude"];
    maxLon = mapgrp["MaximumLongitude"];
    if ((maxLat - minLat) >= 180.0 && (maxLon - minLon) >= 360.0) isGlobal = true;
  }

  Projection *proj = icube->getProjection();
  if(proj == NULL) {
    iString message = "The input cube must be a DEM file, which means it must be projected. ";
    message += "This file is not map projected.";
    throw iException::Message(iException::User, message, _FILEINFO_);
  }

  if(!proj->IsEquatorialCylindrical()) {
    ocube = p.SetOutputCube("TO");
    p.StartProcess(GetStats);

    PvlGroup demRange("Results");
    demRange += PvlKeyword("MinimumRadius", inCubeStats.Minimum(), "meters");
    demRange += PvlKeyword("MaximumRadius", inCubeStats.Maximum(), "meters");
    Application::Log(demRange);

    // Store min/max radii values in new ShapeModelStatistics table
    string shp_name = "ShapeModelStatistics";
    TableField fmin("MinimumRadius",Isis::TableField::Double);
    TableField fmax("MaximumRadius",Isis::TableField::Double);

    TableRecord record;
    record += fmin;
    record += fmax;

    Table table(shp_name,record);
 
    record[0] = Distance(inCubeStats.Minimum(),
                         Distance::Meters).GetKilometers();
    record[1] = Distance(inCubeStats.Maximum(),
                         Distance::Meters).GetKilometers();
    table += record;

//    ocube = new Cube();
//    ocube->open(Filename(ui.GetFilename("TO")).Expanded(), "rw");
    ocube->write(table);
    p.EndProcess();
    return;
  }

  if (proj->LatitudeTypeString() != "Planetocentric") {
    iString message = "The input cube must have Planetocentric latitude type.";
    throw iException::Message(iException::User, message, _FILEINFO_);
  }

  // Determine if the file is global
  bool isPadded = false;
  int insideImage = 0;
  if (!hasExtents) {
    if (proj->LongitudeDomainString() == "360") {
      proj->SetGround(90.0,0.0);
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(90.0,360.0);
      }
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(-90.0,0.0);
      }
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(-90.0,360.0);
      }
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
      }
      if (proj->IsGood() && insideImage == 4) isGlobal = true;
    } else {
      proj->SetGround(90.0,-180.0);
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(90.0,180.0);
      }
      if (proj->IsGood()) {
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(-90.0,-180.0);
      }
      if (proj->IsGood()) { 
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
        proj->SetGround(-90.0,180.0);
      }
      if (proj->IsGood()) { 
        if (proj->WorldX() > 0.0 && proj->WorldX() < ins+1 &&
            proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
          insideImage = insideImage + 1;
        }
      }
      if (proj->IsGood() && insideImage == 4) isGlobal = true;
    }
  }

  if (isGlobal) {
    if (proj->LongitudeDomainString() == "360") {
      if (proj->LongitudeDirectionString() == "PositiveEast") {
        proj->SetGround(90.0,0.0);
      } else {
        proj->SetGround(90.0,360.0);
      }
      if (proj->WorldX() >= 1.0 && proj->WorldY() >= 1.0) isPadded = true;
    } else {
      if (proj->LongitudeDirectionString() == "PositiveEast") {
        proj->SetGround(90.0,-180.0);
      } else {
        proj->SetGround(90.0,180.0);
      }
      if (proj->WorldX() >= 1.0 && proj->WorldY() >= 1.0) isPadded = true;
    }
  }

  // If the file isn't global, then determine if it contains either
  // the south or north pole
  bool hasSPole = false;
  bool hasNPole = false;
  if (!isGlobal) {
    proj->SetGround(90.0,0.0);
    if (proj->IsGood()) {
      if (proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
        hasNPole = true;
        if (proj->WorldY() >= 1.0) isPadded = true;
      }
    }
    proj->SetGround(-90.0,0.0);
    if (proj->IsGood()) {
      if (proj->WorldY() > 0.0 && proj->WorldY() < inl+1) {
        hasSPole = true;
        if (proj->WorldY() <= inl) isPadded = true;
      }
    }
  }

  // Set the padding parameters
  leftPad = 0;
  rightPad = 0;
  topPad = 0;
  bottomPad = 0;

  if (isGlobal && !isPadded) {
    leftPad = 1;
    rightPad = 1;
    topPad = 1;
    bottomPad = 1;
  }
  if (!isGlobal && !isPadded) {
    leftPad = 0;
    rightPad = 0;
    if (hasNPole) {
      topPad = 1;
    } else {
      topPad = 0;
    }
    if (hasSPole) {
      bottomPad = 1;
    } else {
      bottomPad = 0;
    }
  }

  // Compute the output size
  int ns = ins + leftPad + rightPad;
  int nl = inl + topPad + bottomPad;
  int nb = inb;

  double upperLeftCorner = mapgrp["UpperLeftCornerX"];
  upperLeftCorner -= leftPad * proj->Resolution();
  mapgrp.AddKeyword(PvlKeyword("UpperLeftCornerX", upperLeftCorner, "meters"),
                    Pvl::Replace);

  upperLeftCorner = mapgrp["UpperLeftCornerY"];
  upperLeftCorner += topPad * proj->Resolution();
  mapgrp.AddKeyword(PvlKeyword("UpperLeftCornerY", upperLeftCorner, "meters"),
                    Pvl::Replace);


  p.SetOutputCube("TO", ns, nl, nb);
  // Make sure everything is propagated and closed
  p.EndProcess();

  // Now we'll really be processing our input cube
  p.SetInputCube("FROM");

  // We need to create the output file
  ocube = new Cube();
  ocube->open(Filename(ui.GetFilename("TO")).Expanded(), "rw");

  p.StartProcess(DoWrap);

  // Update mapping grp
  ocube->putGroup(mapgrp);

  PvlGroup demRange("Results");
  demRange += PvlKeyword("MinimumRadius", outCubeStats.Minimum(), "meters");
  demRange += PvlKeyword("MaximumRadius", outCubeStats.Maximum(), "meters");
  Application::Log(demRange);

  // Store min/max radii values in new ShapeModelStatistics table
  string shp_name = "ShapeModelStatistics";
  TableField fmin("MinimumRadius",Isis::TableField::Double);
  TableField fmax("MaximumRadius",Isis::TableField::Double);

  TableRecord record;
  record += fmin;
  record += fmax;

  Table table(shp_name,record);

  record[0] = Distance(outCubeStats.Minimum(),
                       Distance::Meters).GetKilometers();
  record[1] = Distance(outCubeStats.Maximum(),
                       Distance::Meters).GetKilometers();
  table += record;

  ocube->write(table);

  p.EndProcess();
  ocube->close();
  delete ocube;
}

void GetStats(Buffer &in, Buffer &out) {
  inCubeStats.AddData(&in[0], in.size());
  for (int i=0; i<in.size(); i++) {
    out[i] = in[i];
  }
}

void DoWrap(Buffer &in) {
  LineManager outMan(*ocube);

  int inputSize = in.size();
  outMan.SetLine(in.Line() + topPad);
  int outputSize = outMan.size();

  Statistics inputLineStats;

  // We only need stats for the first and last lines of the input cube...
  if ((topPad == 1 && in.Line() == 1) || (bottomPad == 1 && in.Line() == inl)) {
    for(int outputIndex = 0; outputIndex < outputSize; outputIndex++) {
      int inputIndex = outputIndex - leftPad;
      if(inputIndex < 0) {
        inputLineStats.AddData(in[inputIndex + inputSize]);
      }
      else if(inputIndex < inputSize) {
        inputLineStats.AddData(in[inputIndex]);
      }
      else {
        inputLineStats.AddData(in[inputIndex - inputSize]);
      }
    }
  }

  // Write top pad?
  if (topPad == 1 && in.Line() == 1) {
    double average = inputLineStats.Average();
    for(int outputIndex = 0; outputIndex < outputSize; outputIndex++) {
      int inputIndex = outputIndex - leftPad;
      if(inputIndex < 0) {
        outMan[outputIndex] = 2.0 * average - in[inputIndex + inputSize];
      }
      else if(inputIndex < inputSize) {
        outMan[outputIndex] = 2.0 * average - in[inputIndex];
      }
      else {
        outMan[outputIndex] = 2.0 * average - in[inputIndex - inputSize];
      }
    }
    outMan.SetLine(1);

    outCubeStats.AddData(&outMan[0], outMan.size());
    ocube->write(outMan);
    outMan.SetLine(2);
  }

  // Copy most data
  for(int outputIndex = 0; outputIndex < outputSize; outputIndex++) {
    int inputIndex = outputIndex - leftPad;
    if(inputIndex < 0) {
      outMan[outputIndex] = in[inputIndex + inputSize];
    }
    else if(inputIndex < inputSize) {
      outMan[outputIndex] = in[inputIndex];
    }
    else {
      outMan[outputIndex] = in[inputIndex - inputSize];
    }
  }

  outCubeStats.AddData(&outMan[0], outMan.size());
  ocube->write(outMan);

  // Write bottom pad?
  if (bottomPad == 1 && in.Line() == inl) {
    double average = inputLineStats.Average();
    for(int outputIndex = 0; outputIndex < outputSize; outputIndex++) {
      int inputIndex = outputIndex - leftPad;
      if(inputIndex < 0) {
        outMan[outputIndex] = 2.0 * average - in[inputIndex + inputSize];
      }
      else if(inputIndex < inputSize) {
        outMan[outputIndex] = 2.0 * average - in[inputIndex];
      }
      else {
        outMan[outputIndex] = 2.0 * average - in[inputIndex - inputSize];
      }
    }
    outMan.SetLine(inl + topPad + bottomPad);
    outCubeStats.AddData(&outMan[0], outMan.size());
    ocube->write(outMan);
  }
}
