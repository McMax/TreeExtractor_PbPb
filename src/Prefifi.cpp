#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cstdio>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "Prefifi.h"
#include "RootWriter.h"
#include "Event.h"
#include "Particle.h"

using namespace std;

void mainanalyze(TTree *particletree, const float beam_momentum, const TString output_filename="Extracted_distributions.root")
{
	cout << "Beta calculated for nucleon mass: " << nucleon_mass << " GeV/c^2 and beam momentum: " << beam_momentum << endl;

	float angle,
		  p1, p2,
		  pt1, pt2,
		  pz_cms1, pz_cms2,
		  E1, E2,
		  E_prot,
		  inv_mass,
		  gbE1, gbE2,
		  theta1, theta2,
		  y1, y2,
		  y_prot_cms,
		  eta1, eta2,
		  angle_j,
		  angle_diff,
		  y_diff,
		  eta_diff;

	bool	positive,
			positive_j;

	int	n[3];
	unsigned int all_particles=0;
	UInt_t	i,j;

	TLorentzVector v1, v2, v;

	unsigned correlations = 0, pos_correlations = 0, neg_correlations = 0, all_correlations = 0, unlike_correlations = 0;

	Event *event = new Event();
	Particle *particleA, *particleB;
	particletree->SetBranchAddress("event",&event);
	Long64_t treeNentries = particletree->GetEntries();
	cout << "Number of events: " << treeNentries << endl;
	Long64_t ev;

	Particles particles;
	Histos histos;
	TFile *root_output_file;

	histos.init();
	particles.init(&histos, beam_momentum);
	particles.newEvent(true);
	root_output_file = new TFile(output_filename,"recreate");

	cout << "Writing events" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		particletree->GetEntry(ev);

		n[Neg] = n[All] = n[Pos] = 0;

		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);

			pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
			p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
			E1 = TMath::Sqrt(pion_mass*pion_mass+p1*p1);
			E_prot = TMath::Sqrt(proton_mass*proton_mass+p1*p1);
			y_prot_cms = 0.5*TMath::Log((E_prot+particleA->GetPz())/(E_prot-particleA->GetPz())) - particles.y_cms;
			v1.SetPxPyPzE(particleA->GetPx(),particleA->GetPy(),particleA->GetPz(),E1);

			//Minimal pT cut
//			if(pt1 < 0.2)
//				continue;

			y1 = 0.5*TMath::Log((E1+particleA->GetPz())/(E1-particleA->GetPz())) - particles.y_cms;
			angle = TMath::ATan2(particleA->GetPy(), particleA->GetPx());

			particles.analyze(particleA,beam_momentum);

			positive = particleA->isPositive();

			if(event->GetNpa() > 1)
			{
				for(j=i+1; j<event->GetNpa(); ++j)
				{
					particleB = event->GetParticle(j);

					pt2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2));
					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));

					//cout << "p1 = " << p1 << " | p2 = " << p2 << endl;

					E2 = TMath::Sqrt(pion_mass*pion_mass+p2*p2);
					E_prot = TMath::Sqrt(proton_mass*proton_mass+p2*p2);
					y_prot_cms = 0.5*TMath::Log((E_prot+particleB->GetPz())/(E_prot-particleB->GetPz())) - particles.y_cms;
					v2.SetPxPyPzE(particleB->GetPx(),particleB->GetPy(),particleB->GetPz(),E2);

					v = v1 + v2;
					inv_mass = v.M();

//					if(inv_mass < 0.285) //GeV dipion (280 MeV) + Coulomb interactions (5 MeV)
//						continue;

					histos.histInvMass->Fill(inv_mass);

					//cout << "E1 = " << E1 << " | E2 = " << E2 << endl;

					gbE1 = particles.calc_gbE(E1);
					gbE2 = particles.calc_gbE(E2);

					//cout << "Beta factor: " << particles.beta << endl;
					//cout << "Gamma factor: " << particles.gamma << endl;
					//cout << "gamma*beta*E1: " << gbE1 << " | gamma*beta*E2: " << gbE2 << endl;

					pz_cms1 = particles.gamma*particleA->GetPz() - gbE1;
					pz_cms2 = particles.gamma*particleB->GetPz() - gbE2;

					//cout << "pz_cms1 = " << pz_cms1 << " | pz_cms2 = " << pz_cms2 << endl;

					y2 = 0.5*TMath::Log((E2+particleB->GetPz())/(E2-particleB->GetPz())) - particles.y_cms;
					
					//Minimal pT cut
//					if(pt2 < 0.2)
//						continue;

					angle_j = TMath::ATan2(particleB->GetPy(), particleB->GetPx());

					//cout << "y1 = " << y1 << " | y2 = " << y2 << endl;

					theta1 = TMath::Abs(TMath::ATan2(pt1,pz_cms1));
					theta2 = TMath::Abs(TMath::ATan2(pt2,pz_cms2));

					//cout << "theta1 = " << theta1 << " | theta2 = " << theta2 << endl;

					eta1 = -TMath::Log(TMath::Tan(0.5*theta1));
					eta2 = -TMath::Log(TMath::Tan(0.5*theta2));

					//cout << "eta1 = " << eta1 << " | eta2 = " << eta2 << endl;
					//cout << "angle1 = " << angle << " | angle2 = " << angle_j << endl;

					positive_j = particleB->isPositive();
					if((angle_diff = TMath::Abs(angle-angle_j)) > TMath::Pi())
						angle_diff = 2*TMath::Pi()-angle_diff;
					y_diff = TMath::Abs(y1-y2);
					eta_diff = TMath::Abs(eta1-eta2);

					histos.histDyDphiAll->Fill(angle_diff, y_diff);
					histos.histDetaDphiAll->Fill(angle_diff, eta_diff);

					++all_correlations;

					if((positive_j == true) && (positive == true))
					{
						++correlations;
						++pos_correlations;

						histos.histDyDphiPos->Fill(angle_diff, y_diff);
						histos.histDetaDphiPos->Fill(angle_diff, eta_diff);
					}
					else if((positive_j == false) && (positive == false))
					{
						++correlations;
						++neg_correlations;
						histos.histDyDphiNeg->Fill(angle_diff, y_diff);
						histos.histDetaDphiNeg->Fill(angle_diff, eta_diff);
					}
					else
					{
						++unlike_correlations;
						histos.histDyDphiUnlike->Fill(angle_diff, y_diff);
						histos.histDetaDphiUnlike->Fill(angle_diff, eta_diff);
					}
				}
			}

			all_particles++;
			n[All]++;

			if(positive)
				n[Pos]++;
			else
				n[Neg]++;

		}

		//cout << "\rEvent " << ev;
		if(!(ev%10))
			cout << "Event " << ev << " / " << treeNentries << endl;

		particles.newEvent();

	}

	//event--;

	cout << endl << "Filling with zeros" << endl;

	cout << "All correlations: " << all_correlations << endl;
	cout << "Like-sign correlations: " << correlations << endl;
	cout << "Positive correlations: " << pos_correlations << endl;
	cout << "Negative correlations: " << neg_correlations << endl;
	cout << "=======================" << endl << "All particles: " << all_particles << ", all events: " << ev << endl;
	cout << "Mean multiplicity: " << (((double)all_particles)/ev) << endl;

	//histos.histCharged->AddBinContent(1,zeros);
	//histos.histChargedNeg->AddBinContent(1,zeros);
	//histos.histChargedPos->AddBinContent(1,zeros);
	//	histos.histCharged->ResetStats();
	//	histos.histChargedNeg->ResetStats();
	//	histos.histChargedPos->ResetStats();
	root_output_file->cd();
	histos.write();
	histos.clear();
	root_output_file->Close();
}
