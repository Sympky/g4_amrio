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
//
/// \file B4bSteppingAction.cc
/// \brief Implementation of the B4bSteppingAction class

#include "B4bSteppingAction.hh"
#include "B4bRunData.hh"
#include "B4DetectorConstruction.hh"
#include "B4bEventAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 B4bSteppingAction::B4bSteppingAction() //( const B4DetectorConstruction* detectorConstruction)
//       : G4UserSteppingAction(),
//         fDetConstruction(detectorConstruction),
//         fEventAction(EventAction)

{
 fDetConstruction = (B4DetectorConstruction*)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 fEventAction = (B4bEventAction*)
               G4RunManager::GetRunManager()->GetUserEventAction();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4bSteppingAction::~B4bSteppingAction()
{
}

//..oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4bSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

    // G4cout<<volume<<"\n";


  // step length
  G4double stepLength = 0.;

  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    // stepLength = step->GetStepLength();
    //G4cout<<"\neasfasfasfasfasgasgfasep = "<<edep;
  // }


  // if(step->GetPreStepPoint()->GetMaterial()->GetName()=="POLYVINILTOLUOL")
  // {fEventAction->AddAbs(edep,stepLength);G4cout<<"\nedep = "<<edep;};
      //
      if ( volume == "scintilator1" ) {
          G4cout<<"scintilator 1 EDEP = "<<edep;
          // fEventAction->AddAbs(edep,stepLength);
       }
      //
       if ( volume == "scintilator2" ) {
          G4cout<<"scintilator 2 EDEP = "<<edep;
         // fEventAction->AddGap(edep,stepLength);
      }
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
