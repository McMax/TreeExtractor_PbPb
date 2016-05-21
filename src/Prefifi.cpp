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

void mainanalyzePb(TTree *particletree, const float energy, const TString output_filename="Extracted_distributions.root")
{
	mainanalyze(particletree, 0, true, energy, "NONE", output_filename);
}

void mainanalyze(TTree *particletree, const int zeros, bool write_to_root, const float energy, string fifivsbpar, const TString output_filename="Extracted_distributions.root")
{
	ofstream prefifi_file("Pre_fifi.txt");
	//ofstream debugfile("Debug.txt");
	ofstream prefifi_b_file[12];
	float 	phi_b[3][12],
		phiSq_b[3][12];
	int	n_b[3][12];
	TString targettype;

	bool with_fifivsbpar = false;
	bool with_prefifi = true;

	if(!(fifivsbpar.compare("NONE")))
	{
		with_fifivsbpar = false;
		with_prefifi = false;
	}
	else if(!(fifivsbpar.compare("PREFIFI")))
	{
		with_fifivsbpar = false;
		with_prefifi = true;
	}
	else if(!(fifivsbpar.compare("EMPTY")))
	{
		with_fifivsbpar = true;
		targettype = "EMPTY";
	}
	else if(!(fifivsbpar.compare("FULL")))
	{
		with_fifivsbpar = true;
		targettype = "FULL";
	}
	else if(!(fifivsbpar.compare("VENUS")))
	{
		with_fifivsbpar = true;
		targettype = "VENUS";
	}
	else if(!(fifivsbpar.compare("VGCALOR")))
	{
		with_fifivsbpar = true;
		targettype = "VGCALOR";
	}
	else
	{
		with_fifivsbpar = false;
		cout << "Target type unknown. Skipping fifivsbpar extraction." << endl;
	}

	if(with_fifivsbpar)
	{
		TString prefifi_b_filename;
		
		cout << "Saving to:" << endl;
		for(int j=0; j<12; j++)
		{
			prefifi_b_filename = "PhiphiVsBpar_";
			prefifi_b_filename += targettype;
			prefifi_b_filename += "/B";
			prefifi_b_filename += j;
			prefifi_b_filename += ".txt";
			cout << prefifi_b_filename << endl;
			prefifi_b_file[j].open(prefifi_b_filename);             
		}
		
	}

	cout << "Beta calculated for nucleon mass: " << nucleon_mass << " GeV/c^2" << endl;

	float 	phi[3],
		phiSq[3],
		angle,
		angle3,

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
		eta_diff,
		tt_distance;
	
	bool	positive,
		positive_j;

	int	n[3];
	UInt_t	i,j;
	UInt_t	pid1, pid2;

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

//	to dE/dx of particles in (deta,dphi) < (0.5,0.5)
	std::set<UInt_t> unique_particles_y;
	std::set<UInt_t> unique_particles_eta;
	std::set<UInt_t> unique_particles_y_025;
	std::set<UInt_t> unique_particles_eta_025;

	if(write_to_root)
	{
		histos.init();
		particles.init(&histos, energy);
		particles.newEvent(true);
		root_output_file = new TFile(output_filename,"recreate");
	}

	cout << "Writing events" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		particletree->GetEntry(ev);
		//debugfile << "Event: " << ev << ": " << event->GetNpa() << " particles" << endl;
		
		phi[Neg] = phi[All] = phi[Pos]= 0.;
		phiSq[Neg] = phiSq[All] = phiSq[Pos] = 0.;
		n[Neg] = n[All] = n[Pos] = 0;

		if(with_fifivsbpar)
		{
			for(j=0; j<12; ++j)
			{
				phi_b[Neg][j] = phi_b[All][j] = phi_b[Pos][j]= 0.;
				phiSq_b[Neg][j] = phiSq_b[All][j] = phiSq_b[Pos][j] = 0.;
				n_b[Neg][j] = n_b[All][j] = n_b[Pos][j] = 0;
			}
		}

		unique_particles_y.clear();
		unique_particles_eta.clear();
		unique_particles_y_025.clear();
		unique_particles_eta_025.clear();

		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);
			pid1 = particleA->GetPid();
			//debugfile << "Particle " << i << ": " << pid1 << " px=" << particleA->GetPx() << endl;

			if((TMath::Abs(particleA->GetBx()) > 4) || (TMath::Abs(particleA->GetBy()) > 2))
				continue;
			pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
			p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
			E1 = TMath::Sqrt(pion_mass*pion_mass+p1*p1);
			E_prot = TMath::Sqrt(proton_mass*proton_mass+p1*p1);
			y_prot_cms = 0.5*TMath::Log((E_prot+particleA->GetPz())/(E_prot-particleA->GetPz())) - particles.y_cms;
			v1.SetPxPyPzE(particleA->GetPx(),particleA->GetPy(),particleA->GetPz(),E1);

//			if(y_prot_cms > (particles.y_cms - 0.5))		//Quick cross-check
//				continue;

			y1 = 0.5*TMath::Log((E1+particleA->GetPz())/(E1-particleA->GetPz())) - particles.y_cms;
			angle = TMath::ATan2(particleA->GetPy(), particleA->GetPx());

			if(write_to_root)
				particles.analyze(particleA,energy);

			//debugfile << particleA->GetPx() << "\t" << particleA->GetPy() << "\t" << particleA->GetPz() << endl;

			//angle = TMath::Sqrt(TMath::Power(particleA->GetPx(),2) + TMath::Power(particleA->GetPy(),2));
			positive = particleA->isPositive();

			//debugfile << (positive ? "1 " : "-1 ") << angle << endl;
			if(!positive)
				angle3 = mk_angle3(angle);

			if(write_to_root && (event->GetNpa() > 1))
			{
				for(j=i+1; j<event->GetNpa(); ++j)
				{
					particleB = event->GetParticle(j);
					pid2 = particleB->GetPid();

					tt_distance = calculate_distance(particleA, particleB);
					histos.histTTAverageDistance->Fill(tt_distance);

				//	if(tt_distance < 0.6)
				//		continue;

					histos.histTTAverageDistanceCut->Fill(tt_distance);

					pt2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2));
					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));

					//cout << "p1 = " << p1 << " | p2 = " << p2 << endl;

					E2 = TMath::Sqrt(pion_mass*pion_mass+p2*p2);
					E_prot = TMath::Sqrt(proton_mass*proton_mass+p2*p2);
					y_prot_cms = 0.5*TMath::Log((E_prot+particleB->GetPz())/(E_prot-particleB->GetPz())) - particles.y_cms;
					v2.SetPxPyPzE(particleB->GetPx(),particleB->GetPy(),particleB->GetPz(),E2);

					v = v1 + v2;
					inv_mass = v.M();
					
					//if(inv_mass < 0.285) //GeV dipion (280 MeV) + Coulomb interactions (5 MeV)
					//	continue;

					histos.histInvMass->Fill(inv_mass);

					/*
					if(y_prot_cms > (particles.y_cms - 0.5))		//Quick cross-check
						continue;
						*/

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

					histos.histDyDphiAll->Fill(angle_diff, (y_diff = TMath::Abs(y1-y2)));
					histos.histDetaDphiAll->Fill(angle_diff, (eta_diff = TMath::Abs(eta1-eta2)));
					
					++all_correlations;

					if((positive_j == true) && (positive == true))
					{
						//debugfile << "E: " << ev << " pid1: " << particleA->GetPid() << " pid2: " << particleB->GetPid() << endl;
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
				//--------------- Pb+Pb 00R - histogram to check unlike-sign correlations near deta-dphi = (0,0)-(0.5,0.5)
						if((angle_diff < 0.5) && (eta_diff < 0.5))
						{
							if(unique_particles_eta.find(pid1) == unique_particles_eta.end())
							{
								histos.histDedx_DetaDphiUnlike_05->Fill(p1,particleA->GetdEdx());
								unique_particles_eta.insert(pid1);
							}
							
							if(unique_particles_eta.find(pid2) == unique_particles_eta.end())
							{
								histos.histDedx_DetaDphiUnlike_05->Fill(p2,particleB->GetdEdx());
								unique_particles_eta.insert(pid2);
							}
						}
						if((angle_diff < 0.5) && (y_diff < 0.5))
						{
							if(unique_particles_y.find(pid1) == unique_particles_y.end())
							{
								histos.histDedx_DyDphiUnlike_05->Fill(p1,particleA->GetdEdx());
								unique_particles_y.insert(pid1);
							}

							if(unique_particles_y.find(pid2) == unique_particles_y.end())
							{
								histos.histDedx_DyDphiUnlike_05->Fill(p2,particleB->GetdEdx());
								unique_particles_y.insert(pid2);
							}
						}
						if((angle_diff < 0.25) && (eta_diff < 0.25))
						{
							if(unique_particles_eta_025.find(pid1) == unique_particles_eta_025.end())
							{
								histos.histDedx_DetaDphiUnlike_025->Fill(p1,particleA->GetdEdx());
								unique_particles_eta_025.insert(pid1);
							}
							
							if(unique_particles_eta_025.find(pid2) == unique_particles_eta_025.end())
							{
								histos.histDedx_DetaDphiUnlike_025->Fill(p2,particleB->GetdEdx());
								unique_particles_eta_025.insert(pid2);
							}
						}
						if((angle_diff < 0.25) && (y_diff < 0.25))
						{
							if(unique_particles_y_025.find(pid1) == unique_particles_y_025.end())
							{
								histos.histDedx_DyDphiUnlike_025->Fill(p1,particleA->GetdEdx());
								unique_particles_y_025.insert(pid1);
							}

							if(unique_particles_y_025.find(pid2) == unique_particles_y_025.end())
							{
								histos.histDedx_DyDphiUnlike_025->Fill(p2,particleB->GetdEdx());
								unique_particles_y_025.insert(pid2);
							}
						}
					}
				}
			}

			n[All]++;
			phi[All] += angle;
			phiSq[All] += angle*angle;

			if(positive)
			{
				n[Pos]++;
				phi[Pos] += angle;
				phiSq[Pos] += angle*angle;
			}
			else
			{
				n[Neg]++;
				phi[Neg] += angle3;
				phiSq[Neg] += angle3*angle3;
			}

			if(with_fifivsbpar)
			{
				for(j=0; j<12; ++j)
				{
					if((TMath::Abs(particleA->GetBx()) < bx[j]) && (TMath::Abs(particleA->GetBy()) < by[j]))
					{
						n_b[All][j]++;
						phi_b[All][j] += angle;
						phiSq_b[All][j] += angle*angle;

						if(positive)
						{
							n_b[Pos][j]++;
							phi_b[Pos][j] += angle;
							phiSq_b[Pos][j] += angle*angle;
						}
						else
						{
							n_b[Neg][j]++;
							phi_b[Neg][j] += angle3;
							phiSq_b[Neg][j] += angle3*angle3;
						}	
					}
				}
			}
		}	

		//debugfile << "----------" << endl;

		if(with_prefifi)
		{
			prefifi_file << 10000 << "\t\t" << n[All] << "\t" << phi[All] << "\t" << phiSq[All] << "\t\t" << 
				n[Neg] << "\t" << phi[Neg] << "\t" << phiSq[Neg] << "\t\t" << 
				n[Pos] << "\t" << phi[Pos] << "\t" << phiSq[Pos] << endl;
		}

		if(with_fifivsbpar)
		{
			for(j=0; j<12; ++j)
			{
				prefifi_b_file[j] << 10000 << "\t\t" << n_b[All][j] << "\t" << phi_b[All][j] << "\t" << phiSq_b[All][j] << "\t\t" << 
					n_b[Neg][j] << "\t" << phi_b[Neg][j] << "\t" << phiSq_b[Neg][j] << "\t\t" << 
					n_b[Pos][j] << "\t" << phi_b[Pos][j] << "\t" << phiSq_b[Pos][j] << endl;
			}
		}

		//cout << "\rEvent " << ev;
		if(!(ev%50))
			cout << "Event " << ev << endl;

		if(write_to_root)
			particles.newEvent();
	}

	//event--;

	cout << endl << "Filling with zeros" << endl;
	int zero_event = 0;
	
	if(with_fifivsbpar)
	{
		while(zero_event+1 <= zeros)
		{
			zero_event++;
			prefifi_file << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;

			for(j=0; j<12; ++j)
				prefifi_b_file[j] << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;
			
			cout << "\rEvent " << zero_event;
		}
		cout << endl << (ev+zero_event) << " lines written to Pre_fifi" << endl;
		prefifi_file.close();
		for(j=0; j<12; ++j)
			prefifi_b_file[j].close();
	}
	else if(with_prefifi)
	{
		while(zero_event+1 <= zeros)
		{
			zero_event++;
			prefifi_file << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;
			cout << "\rEvent " << zero_event;
		}
		cout << endl << (ev+zero_event) << " lines written to Pre_fifi" << endl;
		prefifi_file.close();
	}
	else
		prefifi_file.close();

	cout << "All correlations: " << all_correlations << endl;
	cout << "Like-sign correlations: " << correlations << endl;
	cout << "Positive correlations: " << pos_correlations << endl;
	cout << "Negative correlations: " << neg_correlations << endl;
	//debugfile << "All correlations: " << all_correlations << endl;
	//debugfile << "\nLike-sign correlations: " << correlations << endl;
	//debugfile << "Positive correlations: " << pos_correlations << endl;
	//debugfile << "Negative correlations: " << neg_correlations << endl;

	//debugfile.close();

	if(write_to_root)
	{
		//histos.histCharged->AddBinContent(1,zeros);
		//histos.histChargedNeg->AddBinContent(1,zeros);
		//histos.histChargedPos->AddBinContent(1,zeros);
		histos.histCharged->ResetStats();
		histos.histChargedNeg->ResetStats();
		histos.histChargedPos->ResetStats();
		root_output_file->cd();
		histos.write();
		histos.clear();
		root_output_file->Close();
	}
}

Float_t calculate_distance(Particle* partA, Particle* partB)
{
	static int count;
	static float distance_sum;
	static TVector2 trackA, trackB, mov;

	count = 0;
	distance_sum = 0;

	//VTPC1 start
	if((partA->GetVTPC1_Sx()!=9999) && (partB->GetVTPC1_Sx()!=9999))
	{
		trackA.Set(partA->GetVTPC1_Sx(),partA->GetVTPC1_Sy());
		trackB.Set(partB->GetVTPC1_Sx(),partB->GetVTPC1_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "VTPC1 start. Distance: " << mov.Mod() << endl;
		++count;
	}

	//VTPC1 end
	if((partA->GetVTPC1_Ex()!=9999) && (partB->GetVTPC1_Ex()!=9999))
	{
		trackA.Set(partA->GetVTPC1_Ex(),partA->GetVTPC1_Sy());
		trackB.Set(partB->GetVTPC1_Ex(),partB->GetVTPC1_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "VTPC1 end. Distance: " << mov.Mod() << endl;
		++count;
	}

	//VTPC2 start
	if((partA->GetVTPC2_Sx()!=9999) && (partB->GetVTPC2_Sx()!=9999))
	{
		trackA.Set(partA->GetVTPC2_Sx(),partA->GetVTPC2_Sy());
		trackB.Set(partB->GetVTPC2_Sx(),partB->GetVTPC2_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "VTPC2 start. Distance: " << mov.Mod() << endl;
		++count;
	}

	//VTPC2 end
	if((partA->GetVTPC2_Ex()!=9999) && (partB->GetVTPC2_Ex()!=9999))
	{
		trackA.Set(partA->GetVTPC2_Ex(),partA->GetVTPC2_Sy());
		trackB.Set(partB->GetVTPC2_Ex(),partB->GetVTPC2_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "VTPC2 end. Distance: " << mov.Mod() << endl;
		++count;
	}

	//MTPC start
	if((partA->GetMTPC_Sx()!=9999) && (partB->GetMTPC_Sx()!=9999))
	{
		trackA.Set(partA->GetMTPC_Sx(),partA->GetMTPC_Sy());
		trackB.Set(partB->GetMTPC_Sx(),partB->GetMTPC_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "MTPC start. Distance: " << mov.Mod() << endl;
		++count;
	}

	//MTPC end
	if((partA->GetMTPC_Ex()!=9999) && (partB->GetMTPC_Ex()!=9999))
	{
		trackA.Set(partA->GetMTPC_Ex(),partA->GetMTPC_Sy());
		trackB.Set(partB->GetMTPC_Ex(),partB->GetMTPC_Sy());
		mov = trackB - trackA;
		distance_sum += mov.Mod();
		//cout << "MTPC end. Distance: " << mov.Mod() << endl;
		++count;
	}

	if(count==0)
		return -1;

	distance_sum = distance_sum/count;

	return distance_sum;
}
