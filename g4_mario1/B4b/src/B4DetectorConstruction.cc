#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Tubs.hh"


#include "G4PSTrackLength.hh"
#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
G4ThreadLocal
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr;
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   scintilator1LV(nullptr),
   fgammaPV(nullptr),
   scintilator2LV(nullptr),
   fCheckOverlaps(true)
{
}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
B4DetectorConstruction::~B4DetectorConstruction()
{
}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
//  // Define materials
  DefineMaterials();
//
//  // Define volumes
  return DefineVolumes();
}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
void B4DetectorConstruction::DefineMaterials()
{
//  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
//
//  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  G4double fractionmass;
  G4int natoms, ncomponents;
  G4String name, symbol;



      a = 1.008*g/mole;
      G4Element* elH  = new G4Element(name = "Hidrogen"	,symbol="H" , z= 1., a);

      a = 12.011*g/mole;
      G4Element* elC  = new G4Element(name = "Carbon"  	,symbol="C" , z= 6., a);

  	a = 14.01*g/mole;
      G4Element* elN  = new G4Element(name = "Azot"			,symbol="N" , z= 7., a);

  	a = 16.00*g/mole;
  	G4Element* elO  = new G4Element(name = "Oxigen"  	,symbol="O" , z= 8., a);

  //	a = 55.85*g/mole;
  //	G4Element* elFe = new G4Element(name = "Fier"       ,symbol="Fe", z= 26.,a);

  //	a = 52.00*g/mole;
  //	G4Element* elCr = new G4Element(name = "Crom"		,symbol="Cr",	z= 24.,a);

  //	a = 58.71*g/mole;
  //	G4Element* elNi	= new G4Element(name = "Nichel"		,symbol="Ni",	z= 28.,a);

  //	a = 28.09*g/mole;
  //	G4Element* elSi	= new G4Element(name = "Siliciu"	,symbol="Si",	z= 14.,a);

  //	a = 40.08*g/mole;
  //	G4Element* elCa	= new G4Element(name = "Calciu"		,symbol="Ca",	z= 20.,a);

  //	a = 26.98*g/mole;
  //	G4Element* elAl	= new G4Element(name = "Aluminiu"	,symbol="Al",	z= 13.,a);

  //	a = 23.00*g/mole;
  //	G4Element* elNa	= new G4Element(name = "Sodiu"		,symbol="Na",	z= 11.,a);



    density = 1.290*mg/cm3;
  	G4Material* Aer = new G4Material(name = "Aer" ,density, ncomponents=2);
  	Aer->AddElement(elN, fractionmass = 0.7);
  	Aer->AddElement(elO, fractionmass = 0.3);


    density = 1.032*g/cm3;
    G4Material* TPV = new G4Material(name="POLYVINILTOLUOL"		,density, ncomponents=2);
    TPV->AddElement(elH, natoms=10);
    TPV->AddElement(elC, natoms=9 );




  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
//         // The argon by NIST Manager is a gas with a different density
//
//  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
//
////al pt gamma
  new G4Material("Aluminiu", z=13, a=26.98*g/mole,density= 2.7*g/cm3);
//
//  // Print materials
 G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}
// 2
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
//  // Geometry parameters
//  G4int nofLayers = 10;
 G4double hx = 1.*m;
 G4double hy =  1.*m;
 G4double hz  = 1.*m;

 G4double scint_x = 0.05*m;
 G4double scint_y = 0.05*m;
 G4double scint_z = 0.05*m;


//  auto layerThickness = absoThickness + gapThickness;
//  auto calorThickness = nofLayers * layerThickness;
//  auto worldSizeXY = 1.2 * calorSizeXY;
//  auto worldSizeZ  = 1.2 * calorThickness;
//
//  // Get materials
 auto defaultMaterial = G4Material::GetMaterial("Aer");
// auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
//  auto gapMaterial = G4Material::GetMaterial("liquidArgon");
 auto gammaMaterial = G4Material::GetMaterial("Aluminiu");
 auto TPV = G4Material::GetMaterial("POLYVINILTOLUOL");

//
//  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial || ! gammaMaterial) {
//   G4ExceptionDescription msg;
//   msg << "Cannot retrieve materials already defined.";
//   G4Exception("B4DetectorConstruction::DefineVolumes()",
//     "MyCode0001", FatalException, msg);
//}



//  //
 auto worldS
    = new G4Box("World",           // its name
                 hx*5,hy*5,hz*5); // its size

  auto worldLV
    = new G4LogicalVolume(
                worldS,           // its solid
                defaultMaterial,  // its material
                "World");         // its name

 auto worldPV
   = new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,0,0),  // at (0,0,0)
                worldLV,          // its logical volume
                "World",          // its name
                0,                // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps
//
//  //
//  // Calorimeter




//***********************************************
// -.25 *****************************************
// var scint_xyz cutie scaleaza dupa scint_xyz***
// sursa rad pamant******************************
// tpv scint*************************************
// puta jos 45 cm 15*3 8 grosime*****************
// scint 50 50  3********************************
//***********************************************
////
auto calorimeterSs_plin= new G4Box("Calorimetersus",     // its name
                          5.3*scint_x,5.3*scint_y,0.33*scint_z); // its size
auto calorimeterSs_gol= new G4Box("Calorimetersus",     // its name
                          5.*scint_x,5.*scint_y,0.3*scint_z); // its size
auto  calorimeterSs = new G4SubtractionSolid("box_up", calorimeterSs_plin,
                          calorimeterSs_gol, 0, G4ThreeVector(0.,0., 0.01));
auto calorsLV = new G4LogicalVolume(
                calorimeterSs,     // its solid
                defaultMaterial,  // its material
                "Calorimetersus");   // its name
G4double zpos = -1.*m;
new G4PVPlacement(
                0,                // no rotation
                G4ThreeVector(0,0,zpos),  // at (0,0,0)
                calorsLV,          // its logical volume
                "Calorimetersus",    // its name
                worldLV,          // its mother  volume
                false,            // no boolean operation
                0,                // copy number
                fCheckOverlaps);  // checking overlaps



auto calorimeterSj_plin = new G4Box("Calorimeterjos",     // its name
                          5.3*scint_x,5.3*scint_y,0.33*scint_z); //
auto calorimeterSj_gol= new G4Box("Calorimeterjos",     // its name
                           5.*scint_x,5.*scint_y,0.3*scint_z); // its size
auto  calorimeterSj = new G4SubtractionSolid("box_up", calorimeterSs_plin,
                          calorimeterSs_gol, 0, G4ThreeVector(0.,0., 0.01));
auto calorjLV= new G4LogicalVolume(
               calorimeterSj,     // its solid
               defaultMaterial,  // its material
               "Calorimeterjos");   // its name
G4double zposj = -1.6*m;
new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0,0,zposj),  // at (0,0,0)
               calorjLV,          // its logical volume
               "Calorimeterjos",    // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               0,                // copy number
               fCheckOverlaps);  // checking overlaps





     auto scintilator_1_1 = new G4Box("scintilator1",5.*scint_x,5.*scint_y,0.3*scint_z);
     auto scintilator1LV =  new G4LogicalVolume(
                  scintilator_1_1,     // its solid
                  TPV,  // its material
                  "scintilator1LV");   // its name
     auto scintilator1PV  = new G4PVPlacement(
                                 0,                // no rotation
                                 G4ThreeVector(0,0,0),  // at (0,0,0)
                                 scintilator1LV,          // its logical volume
                                 "scintilator1",          // its name
                                 calorsLV,                // its mother  volume
                                 false,            // no boolean operation
                                 0,                // copy number
                                 fCheckOverlaps);  // checking overlaps
                           //




      auto scintilator_1_2 = new G4Box("scintilator2",5.*scint_x,5.*scint_y,0.3*scint_z);
      auto scintilator2LV =  new G4LogicalVolume(
                   scintilator_1_2,     // its solid
                   TPV,  // its material
                   "scintilator2LV");   // its name
      auto scintilator2PV  = new G4PVPlacement(
                                  0,                // no rotation
                                  G4ThreeVector(0,0,0),  // at (0,0,0)
                                  scintilator2LV,          // its logical volume
                                  "scintilator2",          // its name
                                  calorjLV,                // its mother  volume
                                  false,            // no boolean operation
                                  0,                // copy number
                                  fCheckOverlaps);  // checking overlaps
                                                      //

                           ////

////
//// Layer
//  //
//auto layerS
//  = new G4Box("Layer",           // its name
//               calorSizeXY/2, calorSizeXY/2, layerThickness/2); // its size
//auto layerLV
//  = new G4LogicalVolume(
//               layerS,           // its solid
//               defaultMaterial,  // its material
//               "Layer");         // its nam
//new G4PVReplica(
//               "Layer",          // its name
//               layerLV,          // its logical volume
//               calorLV,          // its mother
//               kZAxis,           // axis of replication
//               nofLayers,        // number of replica
//               layerThickness);  // witdth of replicas
//  //
//  // Absorber
//  //
//  auto absorberS
//    = new G4Box("Abso",            // its name
//                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
//
//  auto absorberLV
//    = new G4LogicalVolume(
//                 absorberS,        // its solid
//                 absorberMaterial, // its material
//                 "Abso");          // its name
//
//  fAbsorberPV
//    = new G4PVPlacement(
//                 0,                // no rotation
//                 G4ThreeVector(0., 0., -gapThickness/2), // its position
//                 absorberLV,       // its logical volume
//                 "Abso",           // its name
//                 layerLV,          // its mother  volume
//                 false,            // no boolean operation
//                 0,                // copy number
//                 fCheckOverlaps);  // checking overlaps
//  //
//  // Gap
//  //
//  auto gapS
//    = new G4Box("Gap",             // its name
//                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size
//
//  auto gapLV
//    = new G4LogicalVolume(
//                 gapS,             // its solid
//                 gapMaterial,      // its material
//                 "Gap");           // its name
//
//  fGapPV
//    = new G4PVPlacement(
//                 0,                // no rotation
//                 G4ThreeVector(0., 0., absoThickness/2), // its position
//                 gapLV,            // its logical volume
//                 "Gap",            // its name
//                 layerLV,          // its mother  volume
//                 false,            // no boolean operation
//                 0,                // copy number
//                 fCheckOverlaps);  // checking overlaps
//
//
//
//  // og gamma
//  //
//  //
  G4double rmin = 0.;
  G4double rmax = 8.*cm;

  G4double phimin = 0.;
  G4double dphi = 360.*deg;
//  G4RotationMatrix* zRot8 = new G4RotationMatrix;
//  G4double zRot = 90.*deg;
  //

//G4RotationMatrix* rotmmp  = new G4RotationMatrix();
//rotmmp.rotateZ(90*deg);
//G4ThreeVector rotatie1p(0, 0, 0);
//G4Transform3D rot1p = G4Transform3D(rotmmp,rotatie1p);

// G4RotationMatrix* rotD3 = new G4RotationMatrix();
// rotD3->rotateY(90.*deg);
G4double lung = 0.225*m;
G4Tubs*gammaS
  = new G4Tubs ("Gamma",
                 rmin, rmax, lung, phimin, dphi);
auto fgammaLV
  = new G4LogicalVolume(
                gammaS,        // its solid
                gammaMaterial, // its material
                "fgammaLV");          // its name
G4double zposg = -0.4*m;
fgammaPV
  = new G4PVPlacement(
               0,                // no rotation
               G4ThreeVector(0., 0.,zposg), // its position
               fgammaLV,       // its logical volume
               "fgammaPV",           // its name
               worldLV,          // its mother  volume
               false,            // no boolean operation
               0,                // copy number
               fCheckOverlaps);  // checking overlaps
//
//
//
//
//
//
//  //
//  // print parameters
//  //
//  G4cout
//    << G4endl
//    << "------------------------------------------------------------" << G4endl
//    << "---> The calorimeter is " << nofLayers << " layers of: [ "
//    << absoThickness/mm << "mm of " << absorberMaterial->GetName()
//    << " + "
//    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
//    << "------------------------------------------------------------" << G4endl;
//
//  //
//  // Visualization attributes
//  //
worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
//
//  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
//  simpleBoxVisAtt->SetVisibility(true);
//  calorLV->SetVisAttributes(simpleBoxVisAtt);
//
//  //
//  // Always return the physical World
//  //
 return worldPV;
}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//void B4DetectorConstruction::ConstructSDandField()
//{
//  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
//
//
//  auto scint1Detector = new G4MultiFunctionalDetector("scintilator1");
//  G4SDManager::GetSDMpointer()->AddNewDetector(scint1Detector);
//
//  G4VPrimitiveScorer* primitive;
//  primitive = new G4PSEnergyDeposit("Edep");
//  scint1Detector->RegisterPrimitive(primitive);
//
//  primitive = new G4PSTrackLength("TrackLength");
//  auto charged = new G4SDChargedFilter("chargedFilter");
//  primitive ->SetFilter(charged);
//  scint1Detector->RegisterPrimitive(primitive);
//
//  SetSensitiveDetector("scintilator1LV",scint1Detector);
//
//  // declare Gap as a MultiFunctionalDetector scorer
//  //
//  auto scint2Detector = new G4MultiFunctionalDetector("scintilator2");
//  G4SDManager::GetSDMpointer()->AddNewDetector(scint2Detector);
//
//  primitive = new G4PSEnergyDeposit("Edep");
//  scint2Detector->RegisterPrimitive(primitive);
//
//  primitive = new G4PSTrackLength("TrackLength");
//  primitive ->SetFilter(charged);
//  scint2Detector->RegisterPrimitive(primitive);
//
//  SetSensitiveDetector("scintilator2LV",scint2Detector);
//  // Create global magnetic field messenger.
//  // Uniform magnetic field is then created automatically if
//  // the field value is not zero.
//  G4ThreeVector fieldValue;
//  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
//  fMagFieldMessenger->SetVerboseLevel(1);
//
//  // Register the field messenger for deleting
//  G4AutoDelete::Register(fMagFieldMessenger);
//}
//
////....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
