//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
// Added by O.Mete
#include "G4ParticleGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 25*mm, env_sizeZ = 2*mm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_W");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 5*env_sizeXY;
  G4double world_sizeZ  = 100*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        env_sizeXY, env_sizeXY, env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Plasma Section
  // 
  // Material definition by using NIST database 
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Galactic");
 
  // Holes /////////////////////////////////////////////
    
    // //////////////////////////////////////////////////////////////////////////////
    
    //  CENTRAL GRID ***************************************************************
    
    // ////////////////////// CENTRE ////////////////////////////////////////////////
  G4double dx = 2.3*mm;
  G4double dy = 2.3*mm;
    
    // Centre
    G4ThreeVector pos1 = G4ThreeVector(0, 0, 1*mm);
    //
    G4ThreeVector pos2 = G4ThreeVector(1*dx, 0, 1*mm);
    G4ThreeVector pos3 = G4ThreeVector(2*dx, 0, 1*mm);
    G4ThreeVector pos4 = G4ThreeVector(3*dx, 0, 1*mm);
    G4ThreeVector pos5 = G4ThreeVector(4*dx, 0, 1*mm);
    G4ThreeVector pos6 = G4ThreeVector(5*dx, 0, 1*mm);
    G4ThreeVector pos7 = G4ThreeVector(6*dx, 0, 1*mm);
    G4ThreeVector pos8 = G4ThreeVector(7*dx, 0, 1*mm);
    G4ThreeVector pos9 = G4ThreeVector(8*dx, 0, 1*mm);
    G4ThreeVector pos10 = G4ThreeVector(9*dx, 0, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n = G4ThreeVector(-1*dx, 0, 1*mm);
    G4ThreeVector pos3n = G4ThreeVector(-2*dx, 0, 1*mm);
    G4ThreeVector pos4n = G4ThreeVector(-3*dx, 0, 1*mm);
    G4ThreeVector pos5n = G4ThreeVector(-4*dx, 0, 1*mm);
    G4ThreeVector pos6n = G4ThreeVector(-5*dx, 0, 1*mm);
    G4ThreeVector pos7n = G4ThreeVector(-6*dx, 0, 1*mm);
    G4ThreeVector pos8n = G4ThreeVector(-7*dx, 0, 1*mm);
    G4ThreeVector pos9n = G4ThreeVector(-8*dx, 0, 1*mm);
    G4ThreeVector pos10n = G4ThreeVector(-9*dx, 0, 1*mm);

  
  G4double shape1_rmin =  0.*mm;
  G4double shape1_rmax =  0.4*mm;
  G4double shape1_hz = 5*mm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Tubs* solidShape1 = new G4Tubs("Shape1", 
                                    shape1_rmin, 
                                    shape1_rmax, 
                                    shape1_hz,
                                    shape1_phimin, 
                                    shape1_phimax);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                            shape1_mat,          //its material
                           "Shape1");           //its name
// /////////////////////////////////////////////////////////////////////////////

  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
// ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 1
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos11 = G4ThreeVector(0, 1*dy, 1*mm);
    //
    G4ThreeVector pos21 = G4ThreeVector(1*dx, 1*dy, 1*mm);
    G4ThreeVector pos31 = G4ThreeVector(2*dx, 1*dy, 1*mm);
    G4ThreeVector pos41 = G4ThreeVector(3*dx, 1*dy, 1*mm);
    G4ThreeVector pos51 = G4ThreeVector(4*dx, 1*dy, 1*mm);
    G4ThreeVector pos61 = G4ThreeVector(5*dx, 1*dy, 1*mm);
    G4ThreeVector pos71 = G4ThreeVector(6*dx, 1*dy, 1*mm);
    G4ThreeVector pos81 = G4ThreeVector(7*dx, 1*dy, 1*mm);
    G4ThreeVector pos91 = G4ThreeVector(8*dx, 1*dy, 1*mm);
    G4ThreeVector pos101 = G4ThreeVector(9*dx, 1*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n1 = G4ThreeVector(-1*dx, 1*dy, 1*mm);
    G4ThreeVector pos3n1 = G4ThreeVector(-2*dx, 1*dy, 1*mm);
    G4ThreeVector pos4n1 = G4ThreeVector(-3*dx, 1*dy, 1*mm);
    G4ThreeVector pos5n1 = G4ThreeVector(-4*dx, 1*dy, 1*mm);
    G4ThreeVector pos6n1 = G4ThreeVector(-5*dx, 1*dy, 1*mm);
    G4ThreeVector pos7n1 = G4ThreeVector(-6*dx, 1*dy, 1*mm);
    G4ThreeVector pos8n1 = G4ThreeVector(-7*dx, 1*dy, 1*mm);
    G4ThreeVector pos9n1 = G4ThreeVector(-8*dx, 1*dy, 1*mm);
    G4ThreeVector pos10n1 = G4ThreeVector(-9*dx, 1*dy, 1*mm);

    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos11,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos21,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos31,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos41,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos51,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos61,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos71,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos81,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos91,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos101,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n1,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 2
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos12 = G4ThreeVector(0, 2*dy, 1*mm);
    //
    G4ThreeVector pos22 = G4ThreeVector(1*dx, 2*dy, 1*mm);
    G4ThreeVector pos32 = G4ThreeVector(2*dx, 2*dy, 1*mm);
    G4ThreeVector pos42 = G4ThreeVector(3*dx, 2*dy, 1*mm);
    G4ThreeVector pos52 = G4ThreeVector(4*dx, 2*dy, 1*mm);
    G4ThreeVector pos62 = G4ThreeVector(5*dx, 2*dy, 1*mm);
    G4ThreeVector pos72 = G4ThreeVector(6*dx, 2*dy, 1*mm);
    G4ThreeVector pos82 = G4ThreeVector(7*dx, 2*dy, 1*mm);
    G4ThreeVector pos92 = G4ThreeVector(8*dx, 2*dy, 1*mm);
    G4ThreeVector pos102 = G4ThreeVector(9*dx, 2*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n2 = G4ThreeVector(-1*dx, 2*dy, 1*mm);
    G4ThreeVector pos3n2 = G4ThreeVector(-2*dx, 2*dy, 1*mm);
    G4ThreeVector pos4n2 = G4ThreeVector(-3*dx, 2*dy, 1*mm);
    G4ThreeVector pos5n2 = G4ThreeVector(-4*dx, 2*dy, 1*mm);
    G4ThreeVector pos6n2 = G4ThreeVector(-5*dx, 2*dy, 1*mm);
    G4ThreeVector pos7n2 = G4ThreeVector(-6*dx, 2*dy, 1*mm);
    G4ThreeVector pos8n2 = G4ThreeVector(-7*dx, 2*dy, 1*mm);
    G4ThreeVector pos9n2 = G4ThreeVector(-8*dx, 2*dy, 1*mm);
    G4ThreeVector pos10n2 = G4ThreeVector(-9*dx, 2*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos12,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos22,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos32,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos42,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos52,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos62,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos72,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos82,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos92,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos102,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n2,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 3
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos13 = G4ThreeVector(0, 3*dy, 1*mm);
    //
    G4ThreeVector pos23 = G4ThreeVector(1*dx, 3*dy, 1*mm);
    G4ThreeVector pos33 = G4ThreeVector(2*dx, 3*dy, 1*mm);
    G4ThreeVector pos43 = G4ThreeVector(3*dx, 3*dy, 1*mm);
    G4ThreeVector pos53 = G4ThreeVector(4*dx, 3*dy, 1*mm);
    G4ThreeVector pos63 = G4ThreeVector(5*dx, 3*dy, 1*mm);
    G4ThreeVector pos73 = G4ThreeVector(6*dx, 3*dy, 1*mm);
    G4ThreeVector pos83 = G4ThreeVector(7*dx, 3*dy, 1*mm);
    G4ThreeVector pos93 = G4ThreeVector(8*dx, 3*dy, 1*mm);
    G4ThreeVector pos103 = G4ThreeVector(9*dx, 3*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n3 = G4ThreeVector(-1*dx, 3*dy, 1*mm);
    G4ThreeVector pos3n3 = G4ThreeVector(-2*dx, 3*dy, 1*mm);
    G4ThreeVector pos4n3 = G4ThreeVector(-3*dx, 3*dy, 1*mm);
    G4ThreeVector pos5n3 = G4ThreeVector(-4*dx, 3*dy, 1*mm);
    G4ThreeVector pos6n3 = G4ThreeVector(-5*dx, 3*dy, 1*mm);
    G4ThreeVector pos7n3 = G4ThreeVector(-6*dx, 3*dy, 1*mm);
    G4ThreeVector pos8n3 = G4ThreeVector(-7*dx, 3*dy, 1*mm);
    G4ThreeVector pos9n3 = G4ThreeVector(-8*dx, 3*dy, 1*mm);
    G4ThreeVector pos10n3 = G4ThreeVector(-9*dx, 3*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos13,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos23,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos33,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos43,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos53,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos63,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos73,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos83,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos93,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos103,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n3,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 4
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos14 = G4ThreeVector(0, 4*dy, 1*mm);
    //
    G4ThreeVector pos24 = G4ThreeVector(1*dx, 4*dy, 1*mm);
    G4ThreeVector pos34 = G4ThreeVector(2*dx, 4*dy, 1*mm);
    G4ThreeVector pos44 = G4ThreeVector(3*dx, 4*dy, 1*mm);
    G4ThreeVector pos54 = G4ThreeVector(4*dx, 4*dy, 1*mm);
    G4ThreeVector pos64 = G4ThreeVector(5*dx, 4*dy, 1*mm);
    G4ThreeVector pos74 = G4ThreeVector(6*dx, 4*dy, 1*mm);
    G4ThreeVector pos84 = G4ThreeVector(7*dx, 4*dy, 1*mm);
    G4ThreeVector pos94 = G4ThreeVector(8*dx, 4*dy, 1*mm);
    G4ThreeVector pos104 = G4ThreeVector(9*dx, 4*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n4 = G4ThreeVector(-1*dx, 4*dy, 1*mm);
    G4ThreeVector pos3n4 = G4ThreeVector(-2*dx, 4*dy, 1*mm);
    G4ThreeVector pos4n4 = G4ThreeVector(-3*dx, 4*dy, 1*mm);
    G4ThreeVector pos5n4 = G4ThreeVector(-4*dx, 4*dy, 1*mm);
    G4ThreeVector pos6n4 = G4ThreeVector(-5*dx, 4*dy, 1*mm);
    G4ThreeVector pos7n4 = G4ThreeVector(-6*dx, 4*dy, 1*mm);
    G4ThreeVector pos8n4 = G4ThreeVector(-7*dx, 4*dy, 1*mm);
    G4ThreeVector pos9n4 = G4ThreeVector(-8*dx, 4*dy, 1*mm);
    G4ThreeVector pos10n4 = G4ThreeVector(-9*dx, 4*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos14,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos24,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos34,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos44,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos54,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos64,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos74,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos84,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos94,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos104,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n4,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 5
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos15 = G4ThreeVector(0, 5*dy, 1*mm);
    //
    G4ThreeVector pos25 = G4ThreeVector(1*dx, 5*dy, 1*mm);
    G4ThreeVector pos35 = G4ThreeVector(2*dx, 5*dy, 1*mm);
    G4ThreeVector pos45 = G4ThreeVector(3*dx, 5*dy, 1*mm);
    G4ThreeVector pos55 = G4ThreeVector(4*dx, 5*dy, 1*mm);
    G4ThreeVector pos65 = G4ThreeVector(5*dx, 5*dy, 1*mm);
    G4ThreeVector pos75 = G4ThreeVector(6*dx, 5*dy, 1*mm);
    G4ThreeVector pos85 = G4ThreeVector(7*dx, 5*dy, 1*mm);
    G4ThreeVector pos95 = G4ThreeVector(8*dx, 5*dy, 1*mm);
    G4ThreeVector pos105 = G4ThreeVector(9*dx, 5*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n5 = G4ThreeVector(-1*dx, 5*dy, 1*mm);
    G4ThreeVector pos3n5 = G4ThreeVector(-2*dx, 5*dy, 1*mm);
    G4ThreeVector pos4n5 = G4ThreeVector(-3*dx, 5*dy, 1*mm);
    G4ThreeVector pos5n5 = G4ThreeVector(-4*dx, 5*dy, 1*mm);
    G4ThreeVector pos6n5 = G4ThreeVector(-5*dx, 5*dy, 1*mm);
    G4ThreeVector pos7n5 = G4ThreeVector(-6*dx, 5*dy, 1*mm);
    G4ThreeVector pos8n5 = G4ThreeVector(-7*dx, 5*dy, 1*mm);
    G4ThreeVector pos9n5 = G4ThreeVector(-8*dx, 5*dy, 1*mm);
    G4ThreeVector pos10n5 = G4ThreeVector(-9*dx, 5*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos15,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos25,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos35,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos45,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos55,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos65,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos75,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos85,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos95,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos105,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n5,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    fScoringVolume = logicShape1;
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 6
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos16 = G4ThreeVector(0, 6*dy, 1*mm);
    //
    G4ThreeVector pos26 = G4ThreeVector(1*dx, 6*dy, 1*mm);
    G4ThreeVector pos36 = G4ThreeVector(2*dx, 6*dy, 1*mm);
    G4ThreeVector pos46 = G4ThreeVector(3*dx, 6*dy, 1*mm);
    G4ThreeVector pos56 = G4ThreeVector(4*dx, 6*dy, 1*mm);
    G4ThreeVector pos66 = G4ThreeVector(5*dx, 6*dy, 1*mm);
    G4ThreeVector pos76 = G4ThreeVector(6*dx, 6*dy, 1*mm);
    G4ThreeVector pos86 = G4ThreeVector(7*dx, 6*dy, 1*mm);
    G4ThreeVector pos96 = G4ThreeVector(8*dx, 6*dy, 1*mm);
    G4ThreeVector pos106 = G4ThreeVector(9*dx, 6*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n6 = G4ThreeVector(-1*dx, 6*dy, 1*mm);
    G4ThreeVector pos3n6 = G4ThreeVector(-2*dx, 6*dy, 1*mm);
    G4ThreeVector pos4n6 = G4ThreeVector(-3*dx, 6*dy, 1*mm);
    G4ThreeVector pos5n6 = G4ThreeVector(-4*dx, 6*dy, 1*mm);
    G4ThreeVector pos6n6 = G4ThreeVector(-5*dx, 6*dy, 1*mm);
    G4ThreeVector pos7n6 = G4ThreeVector(-6*dx, 6*dy, 1*mm);
    G4ThreeVector pos8n6 = G4ThreeVector(-7*dx, 6*dy, 1*mm);
    G4ThreeVector pos9n6 = G4ThreeVector(-8*dx, 6*dy, 1*mm);
    G4ThreeVector pos10n6 = G4ThreeVector(-9*dx, 6*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos16,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos26,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos36,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos46,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos56,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos66,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos76,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos86,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos96,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos106,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n6,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 7
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos17 = G4ThreeVector(0, 7*dy, 1*mm);
    //
    G4ThreeVector pos27 = G4ThreeVector(1*dx, 7*dy, 1*mm);
    G4ThreeVector pos37 = G4ThreeVector(2*dx, 7*dy, 1*mm);
    G4ThreeVector pos47 = G4ThreeVector(3*dx, 7*dy, 1*mm);
    G4ThreeVector pos57 = G4ThreeVector(4*dx, 7*dy, 1*mm);
    G4ThreeVector pos67 = G4ThreeVector(5*dx, 7*dy, 1*mm);
    G4ThreeVector pos77 = G4ThreeVector(6*dx, 7*dy, 1*mm);
    G4ThreeVector pos87 = G4ThreeVector(7*dx, 7*dy, 1*mm);
    G4ThreeVector pos97 = G4ThreeVector(8*dx, 7*dy, 1*mm);
    G4ThreeVector pos107 = G4ThreeVector(9*dx, 7*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n7 = G4ThreeVector(-1*dx, 7*dy, 1*mm);
    G4ThreeVector pos3n7 = G4ThreeVector(-2*dx, 7*dy, 1*mm);
    G4ThreeVector pos4n7 = G4ThreeVector(-3*dx, 7*dy, 1*mm);
    G4ThreeVector pos5n7 = G4ThreeVector(-4*dx, 7*dy, 1*mm);
    G4ThreeVector pos6n7 = G4ThreeVector(-5*dx, 7*dy, 1*mm);
    G4ThreeVector pos7n7 = G4ThreeVector(-6*dx, 7*dy, 1*mm);
    G4ThreeVector pos8n7 = G4ThreeVector(-7*dx, 7*dy, 1*mm);
    G4ThreeVector pos9n7 = G4ThreeVector(-8*dx, 7*dy, 1*mm);
    G4ThreeVector pos10n7 = G4ThreeVector(-9*dx, 7*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos17,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos27,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos37,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos47,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos57,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos67,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos77,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos87,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos97,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos107,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n7,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 8
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos18 = G4ThreeVector(0, 8*dy, 1*mm);
    //
    G4ThreeVector pos28 = G4ThreeVector(1*dx, 8*dy, 1*mm);
    G4ThreeVector pos38 = G4ThreeVector(2*dx, 8*dy, 1*mm);
    G4ThreeVector pos48 = G4ThreeVector(3*dx, 8*dy, 1*mm);
    G4ThreeVector pos58 = G4ThreeVector(4*dx, 8*dy, 1*mm);
    G4ThreeVector pos68 = G4ThreeVector(5*dx, 8*dy, 1*mm);
    G4ThreeVector pos78 = G4ThreeVector(6*dx, 8*dy, 1*mm);
    G4ThreeVector pos88 = G4ThreeVector(7*dx, 8*dy, 1*mm);
    G4ThreeVector pos98 = G4ThreeVector(8*dx, 8*dy, 1*mm);
    G4ThreeVector pos108 = G4ThreeVector(9*dx, 8*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n8 = G4ThreeVector(-1*dx, 8*dy, 1*mm);
    G4ThreeVector pos3n8 = G4ThreeVector(-2*dx, 8*dy, 1*mm);
    G4ThreeVector pos4n8 = G4ThreeVector(-3*dx, 8*dy, 1*mm);
    G4ThreeVector pos5n8 = G4ThreeVector(-4*dx, 8*dy, 1*mm);
    G4ThreeVector pos6n8 = G4ThreeVector(-5*dx, 8*dy, 1*mm);
    G4ThreeVector pos7n8 = G4ThreeVector(-6*dx, 8*dy, 1*mm);
    G4ThreeVector pos8n8 = G4ThreeVector(-7*dx, 8*dy, 1*mm);
    G4ThreeVector pos9n8 = G4ThreeVector(-8*dx, 8*dy, 1*mm);
    G4ThreeVector pos10n8 = G4ThreeVector(-9*dx, 8*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos18,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos28,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos38,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos48,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos58,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos68,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos78,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos88,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos98,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos108,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n8,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 9
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos19 = G4ThreeVector(0, 9*dy, 1*mm);
    //
    G4ThreeVector pos29 = G4ThreeVector(1*dx, 9*dy, 1*mm);
    G4ThreeVector pos39 = G4ThreeVector(2*dx, 9*dy, 1*mm);
    G4ThreeVector pos49 = G4ThreeVector(3*dx, 9*dy, 1*mm);
    G4ThreeVector pos59 = G4ThreeVector(4*dx, 9*dy, 1*mm);
    G4ThreeVector pos69 = G4ThreeVector(5*dx, 9*dy, 1*mm);
    G4ThreeVector pos79 = G4ThreeVector(6*dx, 9*dy, 1*mm);
    G4ThreeVector pos89 = G4ThreeVector(7*dx, 9*dy, 1*mm);
    G4ThreeVector pos99 = G4ThreeVector(8*dx, 9*dy, 1*mm);
    G4ThreeVector pos109 = G4ThreeVector(9*dx, 9*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n9 = G4ThreeVector(-1*dx, 9*dy, 1*mm);
    G4ThreeVector pos3n9 = G4ThreeVector(-2*dx, 9*dy, 1*mm);
    G4ThreeVector pos4n9 = G4ThreeVector(-3*dx, 9*dy, 1*mm);
    G4ThreeVector pos5n9 = G4ThreeVector(-4*dx, 9*dy, 1*mm);
    G4ThreeVector pos6n9 = G4ThreeVector(-5*dx, 9*dy, 1*mm);
    G4ThreeVector pos7n9 = G4ThreeVector(-6*dx, 9*dy, 1*mm);
    G4ThreeVector pos8n9 = G4ThreeVector(-7*dx, 9*dy, 1*mm);
    G4ThreeVector pos9n9 = G4ThreeVector(-8*dx, 9*dy, 1*mm);
    G4ThreeVector pos10n9 = G4ThreeVector(-9*dx, 9*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos19,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos29,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos39,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos49,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos59,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos69,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos79,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos89,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos99,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos109,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n9,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
 
    // //////////////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////
    // NEGATIVE HALF OF THE Y AXIS
    // //////////////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 1
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos11n = G4ThreeVector(0, -1*dy, 1*mm);
    //
    G4ThreeVector pos21n = G4ThreeVector(1*dx, -1*dy, 1*mm);
    G4ThreeVector pos31n = G4ThreeVector(2*dx, -1*dy, 1*mm);
    G4ThreeVector pos41n = G4ThreeVector(3*dx, -1*dy, 1*mm);
    G4ThreeVector pos51n = G4ThreeVector(4*dx, -1*dy, 1*mm);
    G4ThreeVector pos61n = G4ThreeVector(5*dx, -1*dy, 1*mm);
    G4ThreeVector pos71n = G4ThreeVector(6*dx, -1*dy, 1*mm);
    G4ThreeVector pos81n = G4ThreeVector(7*dx, -1*dy, 1*mm);
    G4ThreeVector pos91n = G4ThreeVector(8*dx, -1*dy, 1*mm);
    G4ThreeVector pos101n = G4ThreeVector(9*dx, -1*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n1n = G4ThreeVector(-1*dx, -1*dy, 1*mm);
    G4ThreeVector pos3n1n = G4ThreeVector(-2*dx, -1*dy, 1*mm);
    G4ThreeVector pos4n1n = G4ThreeVector(-3*dx, -1*dy, 1*mm);
    G4ThreeVector pos5n1n = G4ThreeVector(-4*dx, -1*dy, 1*mm);
    G4ThreeVector pos6n1n = G4ThreeVector(-5*dx, -1*dy, 1*mm);
    G4ThreeVector pos7n1n = G4ThreeVector(-6*dx, -1*dy, 1*mm);
    G4ThreeVector pos8n1n = G4ThreeVector(-7*dx, -1*dy, 1*mm);
    G4ThreeVector pos9n1n = G4ThreeVector(-8*dx, -1*dy, 1*mm);
    G4ThreeVector pos10n1n = G4ThreeVector(-9*dx, -1*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos11n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos21n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos31n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos41n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos51n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos61n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos71n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos81n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos91n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos101n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n1n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 2
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos12n = G4ThreeVector(0, -2*dy, 1*mm);
    //
    G4ThreeVector pos22n = G4ThreeVector(1*dx, -2*dy, 1*mm);
    G4ThreeVector pos32n = G4ThreeVector(2*dx, -2*dy, 1*mm);
    G4ThreeVector pos42n = G4ThreeVector(3*dx, -2*dy, 1*mm);
    G4ThreeVector pos52n = G4ThreeVector(4*dx, -2*dy, 1*mm);
    G4ThreeVector pos62n = G4ThreeVector(5*dx, -2*dy, 1*mm);
    G4ThreeVector pos72n = G4ThreeVector(6*dx, -2*dy, 1*mm);
    G4ThreeVector pos82n = G4ThreeVector(7*dx, -2*dy, 1*mm);
    G4ThreeVector pos92n = G4ThreeVector(8*dx, -2*dy, 1*mm);
    G4ThreeVector pos102n = G4ThreeVector(9*dx, -2*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n2n = G4ThreeVector(-1*dx, -2*dy, 1*mm);
    G4ThreeVector pos3n2n = G4ThreeVector(-2*dx, -2*dy, 1*mm);
    G4ThreeVector pos4n2n = G4ThreeVector(-3*dx, -2*dy, 1*mm);
    G4ThreeVector pos5n2n = G4ThreeVector(-4*dx, -2*dy, 1*mm);
    G4ThreeVector pos6n2n = G4ThreeVector(-5*dx, -2*dy, 1*mm);
    G4ThreeVector pos7n2n = G4ThreeVector(-6*dx, -2*dy, 1*mm);
    G4ThreeVector pos8n2n = G4ThreeVector(-7*dx, -2*dy, 1*mm);
    G4ThreeVector pos9n2n = G4ThreeVector(-8*dx, -2*dy, 1*mm);
    G4ThreeVector pos10n2n = G4ThreeVector(-9*dx, -2*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos12n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos22n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos32n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos42n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos52n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos62n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos72n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos82n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos92n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos102n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n2n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 3
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos13n = G4ThreeVector(0, -3*dy, 1*mm);
    //
    G4ThreeVector pos23n = G4ThreeVector(1*dx, -3*dy, 1*mm);
    G4ThreeVector pos33n = G4ThreeVector(2*dx, -3*dy, 1*mm);
    G4ThreeVector pos43n = G4ThreeVector(3*dx, -3*dy, 1*mm);
    G4ThreeVector pos53n = G4ThreeVector(4*dx, -3*dy, 1*mm);
    G4ThreeVector pos63n = G4ThreeVector(5*dx, -3*dy, 1*mm);
    G4ThreeVector pos73n = G4ThreeVector(6*dx, -3*dy, 1*mm);
    G4ThreeVector pos83n = G4ThreeVector(7*dx, -3*dy, 1*mm);
    G4ThreeVector pos93n = G4ThreeVector(8*dx, -3*dy, 1*mm);
    G4ThreeVector pos103n = G4ThreeVector(9*dx, -3*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n3n = G4ThreeVector(-1*dx, -3*dy, 1*mm);
    G4ThreeVector pos3n3n = G4ThreeVector(-2*dx, -3*dy, 1*mm);
    G4ThreeVector pos4n3n = G4ThreeVector(-3*dx, -3*dy, 1*mm);
    G4ThreeVector pos5n3n = G4ThreeVector(-4*dx, -3*dy, 1*mm);
    G4ThreeVector pos6n3n = G4ThreeVector(-5*dx, -3*dy, 1*mm);
    G4ThreeVector pos7n3n = G4ThreeVector(-6*dx, -3*dy, 1*mm);
    G4ThreeVector pos8n3n = G4ThreeVector(-7*dx, -3*dy, 1*mm);
    G4ThreeVector pos9n3n = G4ThreeVector(-8*dx, -3*dy, 1*mm);
    G4ThreeVector pos10n3n = G4ThreeVector(-9*dx, -3*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos13n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos23n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos33n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos43n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos53n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos63n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos73n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos83n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos93n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos103n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n3n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 4
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos14n = G4ThreeVector(0, -4*dy, 1*mm);
    //
    G4ThreeVector pos24n = G4ThreeVector(1*dx, -4*dy, 1*mm);
    G4ThreeVector pos34n = G4ThreeVector(2*dx, -4*dy, 1*mm);
    G4ThreeVector pos44n = G4ThreeVector(3*dx, -4*dy, 1*mm);
    G4ThreeVector pos54n = G4ThreeVector(4*dx, -4*dy, 1*mm);
    G4ThreeVector pos64n = G4ThreeVector(5*dx, -4*dy, 1*mm);
    G4ThreeVector pos74n = G4ThreeVector(6*dx, -4*dy, 1*mm);
    G4ThreeVector pos84n = G4ThreeVector(7*dx, -4*dy, 1*mm);
    G4ThreeVector pos94n = G4ThreeVector(8*dx, -4*dy, 1*mm);
    G4ThreeVector pos104n = G4ThreeVector(9*dx, -4*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n4n = G4ThreeVector(-1*dx, -4*dy, 1*mm);
    G4ThreeVector pos3n4n = G4ThreeVector(-2*dx, -4*dy, 1*mm);
    G4ThreeVector pos4n4n = G4ThreeVector(-3*dx, -4*dy, 1*mm);
    G4ThreeVector pos5n4n = G4ThreeVector(-4*dx, -4*dy, 1*mm);
    G4ThreeVector pos6n4n = G4ThreeVector(-5*dx, -4*dy, 1*mm);
    G4ThreeVector pos7n4n = G4ThreeVector(-6*dx, -4*dy, 1*mm);
    G4ThreeVector pos8n4n = G4ThreeVector(-7*dx, -4*dy, 1*mm);
    G4ThreeVector pos9n4n = G4ThreeVector(-8*dx, -4*dy, 1*mm);
    G4ThreeVector pos10n4n = G4ThreeVector(-9*dx, -4*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos14n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos24n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos34n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos44n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos54n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos64n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos74n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos84n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos94n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos104n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n4n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 5
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos15n = G4ThreeVector(0, -5*dy, 1*mm);
    //
    G4ThreeVector pos25n = G4ThreeVector(1*dx, -5*dy, 1*mm);
    G4ThreeVector pos35n = G4ThreeVector(2*dx, -5*dy, 1*mm);
    G4ThreeVector pos45n = G4ThreeVector(3*dx, -5*dy, 1*mm);
    G4ThreeVector pos55n = G4ThreeVector(4*dx, -5*dy, 1*mm);
    G4ThreeVector pos65n = G4ThreeVector(5*dx, -5*dy, 1*mm);
    G4ThreeVector pos75n = G4ThreeVector(6*dx, -5*dy, 1*mm);
    G4ThreeVector pos85n = G4ThreeVector(7*dx, -5*dy, 1*mm);
    G4ThreeVector pos95n = G4ThreeVector(8*dx, -5*dy, 1*mm);
    G4ThreeVector pos105n = G4ThreeVector(9*dx, -5*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n5n = G4ThreeVector(-1*dx, -5*dy, 1*mm);
    G4ThreeVector pos3n5n = G4ThreeVector(-2*dx, -5*dy, 1*mm);
    G4ThreeVector pos4n5n = G4ThreeVector(-3*dx, -5*dy, 1*mm);
    G4ThreeVector pos5n5n = G4ThreeVector(-4*dx, -5*dy, 1*mm);
    G4ThreeVector pos6n5n = G4ThreeVector(-5*dx, -5*dy, 1*mm);
    G4ThreeVector pos7n5n = G4ThreeVector(-6*dx, -5*dy, 1*mm);
    G4ThreeVector pos8n5n = G4ThreeVector(-7*dx, -5*dy, 1*mm);
    G4ThreeVector pos9n5n = G4ThreeVector(-8*dx, -5*dy, 1*mm);
    G4ThreeVector pos10n5n = G4ThreeVector(-9*dx, -5*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos15n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos25n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos35n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos45n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos55n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos65n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos75n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos85n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos95n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos105n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n5n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    fScoringVolume = logicShape1;
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 6
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos16n = G4ThreeVector(0, -6*dy, 1*mm);
    //
    G4ThreeVector pos26n = G4ThreeVector(1*dx, -6*dy, 1*mm);
    G4ThreeVector pos36n = G4ThreeVector(2*dx, -6*dy, 1*mm);
    G4ThreeVector pos46n = G4ThreeVector(3*dx, -6*dy, 1*mm);
    G4ThreeVector pos56n = G4ThreeVector(4*dx, -6*dy, 1*mm);
    G4ThreeVector pos66n = G4ThreeVector(5*dx, -6*dy, 1*mm);
    G4ThreeVector pos76n = G4ThreeVector(6*dx, -6*dy, 1*mm);
    G4ThreeVector pos86n = G4ThreeVector(7*dx, -6*dy, 1*mm);
    G4ThreeVector pos96n = G4ThreeVector(8*dx, -6*dy, 1*mm);
    G4ThreeVector pos106n = G4ThreeVector(9*dx, -6*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n6n = G4ThreeVector(-1*dx, -6*dy, 1*mm);
    G4ThreeVector pos3n6n = G4ThreeVector(-2*dx, -6*dy, 1*mm);
    G4ThreeVector pos4n6n = G4ThreeVector(-3*dx, -6*dy, 1*mm);
    G4ThreeVector pos5n6n = G4ThreeVector(-4*dx, -6*dy, 1*mm);
    G4ThreeVector pos6n6n = G4ThreeVector(-5*dx, -6*dy, 1*mm);
    G4ThreeVector pos7n6n = G4ThreeVector(-6*dx, -6*dy, 1*mm);
    G4ThreeVector pos8n6n = G4ThreeVector(-7*dx, -6*dy, 1*mm);
    G4ThreeVector pos9n6n = G4ThreeVector(-8*dx, -6*dy, 1*mm);
    G4ThreeVector pos10n6n = G4ThreeVector(-9*dx, -6*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos16n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos26n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos36n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos46n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos56n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos66n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos76n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos86n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos96n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos106n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n6n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   POSITIVE Y 7
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos17n = G4ThreeVector(0, -7*dy, 1*mm);
    //
    G4ThreeVector pos27n = G4ThreeVector(1*dx, -7*dy, 1*mm);
    G4ThreeVector pos37n = G4ThreeVector(2*dx, -7*dy, 1*mm);
    G4ThreeVector pos47n = G4ThreeVector(3*dx, -7*dy, 1*mm);
    G4ThreeVector pos57n = G4ThreeVector(4*dx, -7*dy, 1*mm);
    G4ThreeVector pos67n = G4ThreeVector(5*dx, -7*dy, 1*mm);
    G4ThreeVector pos77n = G4ThreeVector(6*dx, -7*dy, 1*mm);
    G4ThreeVector pos87n = G4ThreeVector(7*dx, -7*dy, 1*mm);
    G4ThreeVector pos97n = G4ThreeVector(8*dx, -7*dy, 1*mm);
    G4ThreeVector pos107n = G4ThreeVector(9*dx, -7*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n7n = G4ThreeVector(-1*dx, -7*dy, 1*mm);
    G4ThreeVector pos3n7n = G4ThreeVector(-2*dx, -7*dy, 1*mm);
    G4ThreeVector pos4n7n = G4ThreeVector(-3*dx, -7*dy, 1*mm);
    G4ThreeVector pos5n7n = G4ThreeVector(-4*dx, -7*dy, 1*mm);
    G4ThreeVector pos6n7n = G4ThreeVector(-5*dx, -7*dy, 1*mm);
    G4ThreeVector pos7n7n = G4ThreeVector(-6*dx, -7*dy, 1*mm);
    G4ThreeVector pos8n7n = G4ThreeVector(-7*dx, -7*dy, 1*mm);
    G4ThreeVector pos9n7n = G4ThreeVector(-8*dx, -7*dy, 1*mm);
    G4ThreeVector pos10n7n = G4ThreeVector(-9*dx, -7*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos17n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos27n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos37n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos47n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos57n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos67n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos77n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos87n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos97n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos107n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n7n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 8
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos18n = G4ThreeVector(0, -8*dy, 1*mm);
    //
    G4ThreeVector pos28n = G4ThreeVector(1*dx, -8*dy, 1*mm);
    G4ThreeVector pos38n = G4ThreeVector(2*dx, -8*dy, 1*mm);
    G4ThreeVector pos48n = G4ThreeVector(3*dx, -8*dy, 1*mm);
    G4ThreeVector pos58n = G4ThreeVector(4*dx, -8*dy, 1*mm);
    G4ThreeVector pos68n = G4ThreeVector(5*dx, -8*dy, 1*mm);
    G4ThreeVector pos78n = G4ThreeVector(6*dx, -8*dy, 1*mm);
    G4ThreeVector pos88n = G4ThreeVector(7*dx, -8*dy, 1*mm);
    G4ThreeVector pos98n = G4ThreeVector(8*dx, -8*dy, 1*mm);
    G4ThreeVector pos108n = G4ThreeVector(9*dx, -8*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n8n = G4ThreeVector(-1*dx, -8*dy, 1*mm);
    G4ThreeVector pos3n8n = G4ThreeVector(-2*dx, -8*dy, 1*mm);
    G4ThreeVector pos4n8n = G4ThreeVector(-3*dx, -8*dy, 1*mm);
    G4ThreeVector pos5n8n = G4ThreeVector(-4*dx, -8*dy, 1*mm);
    G4ThreeVector pos6n8n = G4ThreeVector(-5*dx, -8*dy, 1*mm);
    G4ThreeVector pos7n8n = G4ThreeVector(-6*dx, -8*dy, 1*mm);
    G4ThreeVector pos8n8n = G4ThreeVector(-7*dx, -8*dy, 1*mm);
    G4ThreeVector pos9n8n = G4ThreeVector(-8*dx, -8*dy, 1*mm);
    G4ThreeVector pos10n8n = G4ThreeVector(-9*dx, -8*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos18n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos28n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos38n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos48n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos58n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos68n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos78n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos88n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos98n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos108n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n8n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //   NEGATIVE Y 9
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    // Centre
    G4ThreeVector pos19n = G4ThreeVector(0, -9*dy, 1*mm);
    //
    G4ThreeVector pos29n = G4ThreeVector(1*dx, -9*dy, 1*mm);
    G4ThreeVector pos39n = G4ThreeVector(2*dx, -9*dy, 1*mm);
    G4ThreeVector pos49n = G4ThreeVector(3*dx, -9*dy, 1*mm);
    G4ThreeVector pos59n = G4ThreeVector(4*dx, -9*dy, 1*mm);
    G4ThreeVector pos69n = G4ThreeVector(5*dx, -9*dy, 1*mm);
    G4ThreeVector pos79n = G4ThreeVector(6*dx, -9*dy, 1*mm);
    G4ThreeVector pos89n = G4ThreeVector(7*dx, -9*dy, 1*mm);
    G4ThreeVector pos99n = G4ThreeVector(8*dx, -9*dy, 1*mm);
    G4ThreeVector pos109n = G4ThreeVector(9*dx, -9*dy, 1*mm);
    // NEGATIVE SIDE
    G4ThreeVector pos2n9n = G4ThreeVector(-1*dx, -9*dy, 1*mm);
    G4ThreeVector pos3n9n = G4ThreeVector(-2*dx, -9*dy, 1*mm);
    G4ThreeVector pos4n9n = G4ThreeVector(-3*dx, -9*dy, 1*mm);
    G4ThreeVector pos5n9n = G4ThreeVector(-4*dx, -9*dy, 1*mm);
    G4ThreeVector pos6n9n = G4ThreeVector(-5*dx, -9*dy, 1*mm);
    G4ThreeVector pos7n9n = G4ThreeVector(-6*dx, -9*dy, 1*mm);
    G4ThreeVector pos8n9n = G4ThreeVector(-7*dx, -9*dy, 1*mm);
    G4ThreeVector pos9n9n = G4ThreeVector(-8*dx, -9*dy, 1*mm);
    G4ThreeVector pos10n9n = G4ThreeVector(-9*dx, -9*dy, 1*mm);
    
    // /////////////////////////////////////////////////////////////////////////////
    
    new G4PVPlacement(0,                       //no rotation
                      pos19n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    // ////////////////////// CENTRE ////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos29n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos39n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos49n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos59n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos69n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos79n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos89n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos99n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos109n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    //////////////////////////////////////////////////////////////////////
    // NEGATIVE SIDE /////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    new G4PVPlacement(0,                       //no rotation
                      pos2n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos3n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos4n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos5n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos6n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos7n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos8n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos9n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    new G4PVPlacement(0,                       //no rotation
                      pos10n9n,                    //at position
                      logicShape1,             //its logical volume
                      "Shape1",                //its name
                      logicEnv,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
    

    
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
