#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <TChain.h>
#include <TH1F.h>
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TProof.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "signal.h"
#include <iostream>
#include "TString.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "TROOT.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TGraph.h"
#include "CoincidenceMultiplePmt.hpp"

using namespace std;
std::vector<UCN> Coincidence::SetCoincidence()
{

    int photon_cut = 8.;

    int photon_count = 0.;

    int pmt_hit = 0.;

    UCN event;

    auto itr1_first = std::find_if(v1.begin(), v1.end(), [](const double &x)
                                   { return x > 0; });

    auto itr2_first = std::find_if(v2.begin(), v2.end(), [](const double &x)
                                   { return x > 0; });

    double coincidence_window = 100.;
    double counting_window = 4000.;

    for (auto itr1 = itr1_first; itr1 != v1.end();
         itr1++)
    {
        std::vector<int>::iterator itr2;

        double t1 = *itr1;

        double counting_window1 = t1 + counting_window * 1e-9;
        if (counting_window1 > v1.back())
        {
            break;
        }

        int index_v1_first = 0;
        int index_v1_upper = upper_bound(v1.begin(), v1.end(), t1) - v1.begin();
        int index_v1_lower = lower_bound(v1.begin(), v1.end(), t1) - v1.begin();

        index_v1_first = index_v1_lower;

        int index_v2_first = lower_bound(v2.begin(), v2.end(), t1) - v2.begin();       // look for the hit in second PMT greater than or equal to the one in PMT 1 in time
        int index_v2_first_upper = upper_bound(v2.begin(), v2.end(), t1) - v2.begin(); // look for the hit in second PMT lower than or equal to the one in PMT 2 in time

        double t2 = v2[index_v2_first];
        double t3 = v2[index_v2_first_upper];

        double t4 = 0;
        if (abs(t2 - t1) < abs(t3 - t1))
        {
            t4 = t2;
            index_v2_first = index_v2_first;
        }
        else
        {
            t4 = t3;
            index_v2_first = index_v2_first_upper;
        }

        if (index_v1_first > v1.size() - 1)
        {
            break;
        }

        if (index_v2_first > v2.size() - 1)
        {
            break;
        }

        double t1_next = v1[index_v2_first + 1];
        double t1_next1 = v1[index_v2_first + 2];

        double t2_next = v2[index_v2_first + 1];
        double t2_next1 = v2[index_v2_first + 2];

        std::cout << std::setprecision(19);

        double counting_window2;

        if (t1 > t4)
        {
            pmt_hit = 2;
        }
        else
        {

            pmt_hit = 1;
        }

        counting_window1 = t1 + counting_window * 1e-9; // sec
        counting_window2 = t4 + counting_window * 1e-9; // sec

        if (counting_window2 > v2.back())
        {
            break;
        }

        if (counting_window1 > v1.back())
        {
            break;
        }

        int photon_counter = 0;

        int Na = 0;
        int Nb = 0;

        int index_v1_last = upper_bound(v1.begin(), v1.end(), counting_window1) - v1.begin();
        int index_v2_last = upper_bound(v2.begin(), v2.end(), counting_window2) - v2.begin();

        if (index_v1_last > v1.size() - 1)
        {
            break;
        }

        if (index_v2_last > v2.size() - 1)
        {
            break;
        }

        double diff = (t4 - t1) / 1e-9;

        int NA = 0;
        int NB = 0;

        if (abs(diff) < coincidence_window) // Check for coincidence also check the condition for the event right next to the one found
        {
            double event_length;
            event_length = double(v2.at(index_v2_last) - v1.at(index_v1_first)) * 1e+9;

            NA = index_v1_last - index_v1_first + 1;
            NB = index_v2_last - index_v2_first + 1;

            double weight = event_length * 1e-9;

            photon_count = NA + NB;

            event.StartTime = -999;
            event.EventLength = -999;
            event.ClusterSize = -999;

            if (NA + NB > photon_cut && t1 > 0 && t4 > 0)
            {
                if (t1 < t4)
                {
                    event.StartTime = t1;
                    event.EventLength = event_length;
                    event.ClusterSize = NA + NB;
                }
                else
                {
                    event.StartTime = t4;
                    event.EventLength = event_length;
                    event.ClusterSize = NA + NB;
                }
            }

            CoincidenceVector.push_back(event);

            int advance_index_v1 = index_v1_last - index_v1_first + 1;
            int advance_index_v2 = index_v2_last - index_v2_first + 1;

            if (advance_index_v1 + index_v1_last + 1 > v1.size() || advance_index_v2 + index_v2_last + 1 > v2.size())
            {
                break;
            }

            else
            {

                std::advance(itr1, advance_index_v1);
                std::advance(itr2, advance_index_v2);
            }
        }
        else if (abs(diff) > coincidence_window)
        {
            continue;
        }
    }
    return CoincidenceVector;
}

std::vector<std::pair<double, int>> ElectronicsDeadTime::ElectronicsDeadtimeCorrector()
{
    std::vector<std::pair<double, int>> output;

    sort(pmt_singles.begin(), pmt_singles.end());
    double t0 = pmt_singles[0].first;
    output.push_back(make_pair(t0, pmt_singles[0].second));

    for (int i = 1.; i < pmt_singles.size(); i++)
    {
        if ((pmt_singles[i].first - t0) > 20. * 1e-9)

        {
            t0 = pmt_singles[i].first;
            output.push_back(make_pair(t0, pmt_singles[i].second));
        }
    }

    return output;
}

/*Implementing Lara's Coincidence Algorithm*/
std::vector<UCN> Coincidence::CoincidenceVectorPair()
{
    sort(vector_pair.begin(), vector_pair.end());
    std::vector<UCN> CoincidenceVector;
    double coincidence_window = 100. * 1e-9;
    double tele_window = 1000. * 1e-9;
    double prompt_window = 1000. * 1e-9;
    double Npe = 8.;

    auto ti = vector_pair.begin();

    vector<double> karr;
    vector<double> scalarr;
    vector<double> tauarr;
    vector<double> numphotons;
    vector<double> starttimes;
    vector<double> endtimes;

    UCN event;

    double rate = 0;

    double prevfreephoton = 0;
    double freephotoncount = 0;
    double Npileup;
    srand(time(NULL));
    while (ti != vector_pair.end())
    {
        int nt = 1;       // number of photons found within telescoping window
        int np = 1;       // number of photons found within prompt window
        int state = 0;    // state 0 = no coincidence, state 1 = coincidence found
        int coincevt = 0; // if coincevt = 1, there's enough PEs for an event
        auto tj = ti;
        tj++; // tj starts as ti+1
              // photons->Fill((*tj).realtime-(*ti).realtime);
        while (tj != vector_pair.end())
        {
            if (state == 0)
            {
                if ((*tj).first - (*ti).first >
                    coincidence_window)
                { // break if no PEs found within threshold
                    break;
                }
                else if ((*tj).first - (*ti).first < coincidence_window)
                {
                    nt = nt + 1;
                    np = np + 1;
                    // photons->Fill((*tj).realtime-(*ti).realtime);
                    if ((*tj).second !=
                        (*ti).second)
                    { // state = 1 only if a PE is found in the other
                      // PMT within threshold
                        state = 1;
                        // double ku = (*tj).realtime-(*ti).realtime;
                        // karr.push_back(ku);
                    }
                }
            }
            else if (state == 1)
            {
                auto tj2 = tj;
                --tj2;
                if ((*tj).first - (*tj2).first >
                    tele_window)
                { // continue only if the time between current and last
                  // event is less than telescoping window
                    if (np >= Npe)
                    {
                        coincevt = 1;
                        break;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    if ((*tj).first - (*ti).first < prompt_window)
                    {
                        np = np + 1;
                    }
                    nt = nt + 1;
                    // photons->Fill((*tj).realtime-(*ti).realtime);
                }
            }
            tj++;
        }
        if (coincevt == 1)
        { // if nt reaches threshold, add coincidence event and
          // restart algorithm at current tj
            auto tj2 = tj;
            --tj2;
            double length = ((*tj2).first - (*ti).first);

            /*Fill your coincidence vector here*/

            double randNum = (float)rand() / (float)RAND_MAX;

            // Calculate pileup corrected UCN
            Npileup = nt - length / rate + randNum;
            /*define a UCN event here*/
            event.ClusterSize = nt;
            event.EventLength = length * 1e+9 + 1000.;
            event.Na = 0;
            event.Nb = 0;
            event.StartTime = (*ti).first;
            event.Npile = Npileup;
            CoincidenceVector.push_back(event);
            ti = tj;
        }
        else if (coincevt ==
                 0)
        { // if nt doesn't reach threshold, restart algorithm at ti+1

            // if the photon is not part of a UCN event, use it to calculate the rate
            // of unassociated photons
            if (freephotoncount == 0)
            {
                rate = (*ti).first - prevfreephoton;
            }
            else if (freephotoncount < 100)
            {
                rate = ((freephotoncount - 1) / freephotoncount) * rate +
                       (1 / freephotoncount) * ((*ti).first - prevfreephoton);
            }
            else
            {
                rate = (0.99) * rate + (0.01) * ((*ti).first - prevfreephoton);
            }
            prevfreephoton = (*ti).first;
            freephotoncount = freephotoncount + 1;
            ti++;
        }
    }
    return CoincidenceVector;
}

/****Implementing Frank's Method for Coincidence****/

/*std::vector<UCN> Coincidence::CoincidenceVectorPair()
{

    // int photon_cut = 8.;
    int photon_count = 0.;
    int photon_threshold = 8.0;
    int photon_upper = 70.0;
    double rate = 0;
    double OrphanPhoton = 0;
    double OrphanPhotonTime = 0;
    std::vector<UCN> CoincidenceVector;

    sort(vector_pair.begin(), vector_pair.end());

    UCN event;

    double coincidence_window = 100.;  // in ns
    double telescoping_window = 1000.; // in ns

    event.StartTime = -999;
    event.EventLength = -999;
    event.ClusterSize = -999;
    int initial_pmt = 0;
    event.Na = 0;
    event.Nb = 0;
    event.Npile = 0;

    for (int i = 0; i < vector_pair.size(); i++)
    {

        double event_length = 0.;
        double event_time = 0.;

        int NpE = 1.;
        int k = 0.;
        int NpA = 1.;
        int NpB = 1.;
        int photon_count = 1.;
        double PileCorrectedPhotons = 0;

        double T_diff = (vector_pair[i + 1].first - vector_pair[i].first) * 1e+9;
        int D_diff = vector_pair[i + 1].second - vector_pair[i].second;

        if (T_diff < coincidence_window) // There is a Coincidence

        {

            // std::cout << D_diff << std::endl;
            k = i;
            //  std::cout << "GotHere" << '\t' << i << '\t' << vector_pair.size() << std::endl;
            if (D_diff != 0)
            {
                while ((vector_pair[k + 1].first - vector_pair[k].first) * 1e+9 < telescoping_window)
                {

                    photon_count++;

                    k++;
                    if (vector_pair[k].second == 1)
                    {
                        NpA++;
                    }
                    else if (vector_pair[k].second == 2)
                    {
                        NpB++;
                    }
                    if (k > vector_pair.size())
                    {
                        break;
                    }
                }
            }

            if (photon_count >= photon_threshold)
            {
                if (i + photon_count > vector_pair.size())
                {
                    continue;
                }

                event_length = (vector_pair[i + photon_count].first - vector_pair[i].first) * 1e+9;
                event_time = vector_pair[i].first;

                NpE = photon_count;
                PileCorrectedPhotons = NpE - (event_length * 1e-9) / rate + (float)rand() / (float)RAND_MAX;
                // Photons" << '\t' << NpE << '\t' << PileCorrectedPhotons << std::endl;
                event.StartTime = event_time;
                event.EventLength = event_length + 1000.;
                event.ClusterSize = PileCorrectedPhotons;
                event.Na = NpA;
                event.Nb = NpB;
                event.Npile = PileCorrectedPhotons;
                CoincidenceVector.push_back(event);
                i += (photon_count);
                if (i > vector_pair.size())
                {
                    continue;
                }
            }
        }
        else if (photon_count < photon_threshold)
        {
            if (OrphanPhoton == 0)
            {
                rate = vector_pair[i].first - OrphanPhotonTime;
            }
            else if (OrphanPhoton < 100.)
            {
                rate = ((OrphanPhoton - 1) / OrphanPhoton) * rate +
                       (1 / OrphanPhoton) * (vector_pair[i].first - OrphanPhoton);
            }
            else
            {
                rate = (0.99) * rate + (0.01) * (vector_pair[i].first - OrphanPhotonTime);
            }

            OrphanPhotonTime = vector_pair[i].first;
            OrphanPhotonTime++;
        }
    }
    return CoincidenceVector;
}*/

void Coincidence::PrintElements()
{

    for (int i = 0; i < CoincidenceVector.size(); i++)
    {
        std::cout << CoincidenceVector[i].ClusterSize << std::endl;
        std::cout << CoincidenceVector[i].StartTime << std::endl;
        std::cout << CoincidenceVector[i].EventLength << std::endl;
    }
}

void TTreeWriter::TTreeOutput()
{

    std::string TreeName = "UCN" + std::to_string(det1) + std::to_string(det2);
    TTree UCN(TreeName.c_str(), "UCN event tree");
    double_t length, npe, start, NA, NB, NPILE;
    length = 0;
    npe = 0;
    start = 0;
    NA = 0;
    NB = 0;

    UCN.Branch("start", &start, "start/D");
    UCN.Branch("npe", &npe, "npe/D");
    UCN.Branch("length", &length, "length/D");
    UCN.Branch("NA", &NA, "NA/D");
    UCN.Branch("NB", &NB, "NB/D");
    UCN.Branch("NPILE", &NPILE, "NPILE/D");

    for (int i = 0; i < vector.size(); i++)
    {
        start = vector[i].StartTime;
        npe = vector[i].ClusterSize;
        length = vector[i].EventLength;
        NA = vector[i].Na;
        NB = vector[i].Nb;
        NPILE = vector[i].Npile;
        UCN.Fill();
    }

    UCN.Write();
}

TH1F HistoGramGenerator::HistOutPut()
{

    TH1F h1("h1", "", 40000, 0, 4000);

    for (int i = 0; i < vector.size(); i++)
    {
        h1.Fill(vector[i].StartTime);
    }
    return h1;
}

TH1F HistoGramGenerator::HistOutDead()
{

    TH1F h_dead("h_dead", "", 40000, 0, 4000);

    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i].ClusterSize >= 8.)
        {
            h_dead.Fill(vector[i].StartTime);
        }
    }
    return h_dead;
}

TH1F HistoGramGenerator::HistOutPutDeadWeighted()
{
    TH1F h2_dead("h2_dead", "", 40000, 0, 4000);

    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i].ClusterSize >= 8.)
        {
            h2_dead.Fill(vector[i].StartTime, vector[i].EventLength * 1e-9);
        }
    }
    return h2_dead;
}

TH1F HistoGramGenerator::HistOutPutWeighted()
{
    TH1F h2("h2", "", 40000, 0, 4000);

    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i].Npile >= 8.)
        {
            h2.Fill(vector[i].StartTime, vector[i].EventLength * 1e-9);
        }
    }
    return h2;
}

TH1F HistoGramGenerator::HistOutPutEventlength()
{

    TH1F h3("h3", "", 23000.0, 0.0, 23000.);

    for (int i = 0; i < vector.size(); i++)
    {
        h3.Fill(vector[i].EventLength);
    }
    return h3;
}

TH1F HistoGramGenerator::HistOutPhoton()
{

    TH1F h4("h4", "", 300.0, 0.0, 300.);

    for (int i = 0; i < vector.size(); i++)
    {
        h4.Fill(vector[i].ClusterSize);
    }
    return h4;
}

TH1F DeadTimeCorrection::DeadTimeCorrectedHistogram()
{
    TH1F h_corrected("h_corrected", "", 40000, 0, 4000);

    for (int bin = 1; bin <= h_finer.GetNbinsX(); bin++)
    {

        double counts = h_finer.GetBinContent(bin);
        double weighted_counts = h_weighted.GetBinContent(bin);
        double correction_factor = 0.1 / (0.1 - weighted_counts);
        double final_factor = 0.;
        /*  if (correction_factor < 1.0 || correction_factor > 1.1)
          {
              final_factor = 1.0;
          }

          else
          {
              final_factor = correction_factor;
          }*/

        double corrected_counts = counts * correction_factor;
        h_corrected.SetBinContent(bin, corrected_counts);
    }
    return h_corrected;
}

TH1F ElectronicsDeadTime::TimeDifference()
{

    sort(pmt_singles.begin(), pmt_singles.end());

    TH1F T_Diff("T_Diff", "", 500., 0.0, 500.);

    double t0 = pmt_singles[0].first;

    for (int i = 1.; i < pmt_singles.size(); i++)
    {
        double TDiff = (pmt_singles[i].first - t0) * 1e+9;
        T_Diff.Fill(TDiff);
        t0 = pmt_singles[i].first;
        // std::cout << TDiff << std::endl;
    }
    return T_Diff;
}

TH1F PileUpCorrection::PileUpCorrectedHistogram()
{

    TH1F PileUpCorrectedPMT("PileCorrectedPMT", "", 40000, 0, 4000);

    for (int i = 0; i < vector.size(); i++)
    {
        if (vector[i].Npile >= 8.)
        {
            PileUpCorrectedPMT.Fill(vector[i].StartTime);
        }
    }
    return PileUpCorrectedPMT;
}

/*std::vector<UCN> Coincidence::SetCoincidence()
{

    for (int i = 0; i < v1.size(); i++)
    {
        UCN event;
        event.ClusterSize = v1[i] + v2[i];
        event.EventLength = v1[i] + v2[i];
        event.StartTime = v1[i] + v2[i];
        CoincidenceVector.push_back(event);
    }
    return CoincidenceVector;
}*/

//
int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << "File Name" << std::endl;
        return 1;
    }

    string n1 = argv[1];

    // int file_number = n;

    auto PMT1 = TH1F("PMT1", "PMT1", 4000, 0, 4000);
    auto PMT2 = TH1F("PMT2", "PMT2", 4000, 0, 4000);
    auto PMT3 = TH1F("PMT3", "PMT3", 4000, 0, 4000);
    auto PMT4 = TH1F("PMT4", "PMT4", 4000, 0, 4000);
    auto PMT5 = TH1F("PMT5", "PMT5", 4000, 0, 4000);
    auto PMT6 = TH1F("PMT6", "PMT6", 4000, 0, 4000);
    auto PMT7 = TH1F("PMT7", "PMT7", 4000, 0, 4000);
    auto PMT8 = TH1F("PMT8", "PMT8", 4000, 0, 4000);
    auto h_GV = TH1F("h_GV", "h_GV", 4000, 0, 4000);
    auto h_RH = TH1F("h_RH", "h_RH", 4000, 0, 4000);
    auto h_AC1 = TH1F("h_AC1", "h_AC1", 4000, 0, 4000);
    auto h_AC2 = TH1F("h_AC2", "h_AC2", 4000, 0, 4000);
    auto h_AC3 = TH1F("h_AC3", "h_AC3", 4000, 0, 4000);

    // Open the file containing the tree.

    string file_directory = "/Users/msingh/Desktop/UCNTau/data/replayed_data/";

    // string n1 = std::to_string(n);
    string file_type = "processed_output_";
    string file_suffix = ".root";

    string file_input = file_directory + file_type + n1 + file_suffix;

    string file_output_directory = "/Users/msingh/Desktop/UCNTau/data/coincidence_files/";
    string file_output_type = "processed_output_coincidences_PE_8_offset_check";

    string file_output = file_output_directory + file_output_type + n1 + file_suffix;
    // TFile *myFile = TFile::Open("/Users/msingh/Desktop/UCNTau/data/replayed_data/processed_output_30915.root");

    std::cout << "File Loaded :: " << n1 << std::endl;
    TFile *myFile = TFile::Open(file_input.c_str());

    bool tree1 = myFile->GetListOfKeys()->Contains("tmcs_0");
    bool tree2 = myFile->GetListOfKeys()->Contains("tmcs_1");
    bool tree3 = myFile->GetListOfKeys()->Contains("tmcs_2");

    if (!myFile || myFile->IsZombie())
    {
        std::cout << "Error opening file" << std::endl;
        exit(-1);
    }

    if (tree1 == true && tree2 == true && tree3 == true)
    {
        TTreeReader mcs_0("tmcs_0", myFile);

        TTreeReader mcs_1("tmcs_1", myFile);

        TTreeReader mcs_2("tmcs_2", myFile);

        TTreeReaderValue<double> time(mcs_0, "realtime");
        TTreeReaderValue<int> detnum(mcs_0, "channel");
        TTreeReaderValue<int> tag(mcs_0, "tag");

        TTreeReaderValue<double> time1(mcs_1, "realtime");
        TTreeReaderValue<int> detnum1(mcs_1, "channel");
        TTreeReaderValue<int> tag1(mcs_1, "tag");

        TTreeReaderValue<double> time2(mcs_2, "realtime");
        TTreeReaderValue<int> detnum2(mcs_2, "channel");
        TTreeReaderValue<int> tag2(mcs_2, "tag");

        double timestart1, timestart;

        /* if (!CheckValue(detnum))
             return false;
         if (!CheckValue(time))
             return false;
         if (!CheckValue(detnum1))
             return false;
         if (!CheckValue(time1))
             return false;
         if (!CheckValue(detnum2))
             return false;
         if (!CheckValue(time2))
             return false;*/

        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        int count5 = 0;
        int count6 = 0;
        int count7 = 0;
        int count8 = 0;

        double timing1 = 0;
        double timing2 = 0;
        double timing3 = 0;
        double timing4 = 0;
        double timing5 = 0;
        double timing6 = 0;
        double timing7 = 0;
        double timing8 = 0;
        double timing9 = 0;
        double timing10 = 0;
        double timing11 = 0;
        double timing12 = 0;

        int t1, t2, t3, t4, t5, t6, t7, t8;

        std::vector<double> v1;
        std::vector<double> v2;
        std::vector<double> v3;
        std::vector<double> v4;

        std::vector<double> v5;
        std::vector<double> v6;
        std::vector<double> v7;
        std::vector<double> v8;

        std::vector<std::pair<double, int>> vector_pair1;
        std::vector<std::pair<double, int>> vector_pair2;
        std::vector<std::pair<double, int>> vector_pair3;
        std::vector<std::pair<double, int>> vector_pair4;
        std::vector<std::pair<double, int>> vector_pair23;

        // std::vector<std::vector<std::pair<double, int>>> vec_pairs;

        while (mcs_0.Next())
        {

            if (*tag == 4)
            {
                timestart = *time;
                break;
            }
        }

        while (mcs_1.Next())
        {

            if (*tag1 == 4)
            {
                timestart1 = *time1;
                break;
            }
        }

        while (mcs_0.Next())
        {

            if (*detnum == 1)
            {
                timing1 = *time - timestart;
                if (timing1 > 0)
                {
                    v1.push_back(timing1);
                    PMT1.Fill(timing1);
                    vector_pair1.push_back(make_pair(timing1, 1.));
                }
            }
            else if (*detnum == 2)
            {
                timing2 = *time - timestart;
                if (timing2 > 0)
                {
                    v2.push_back(timing2);
                    PMT2.Fill(timing2);
                    vector_pair1.push_back(make_pair(timing2, 2.));
                    vector_pair23.push_back(make_pair(timing2, 2.));
                }
            }
            else if (*detnum == 3)
            {
                timing3 = *time - timestart;
                if (timing3 > 0)
                {
                    v3.push_back(timing3);
                    PMT3.Fill(timing3);
                    vector_pair2.push_back(make_pair(timing3, 3.));
                    vector_pair23.push_back(make_pair(timing3, 3.));
                }
            }
            else if (*detnum == 4)
            {
                timing4 = *time - timestart;
                if (timing4 > 0)
                {
                    v4.push_back(timing4);
                    PMT4.Fill(timing4);
                    vector_pair2.push_back(make_pair(timing4, 4.));
                }
            }
            else if (*detnum == 5)
            {
                timing10 = *time;
                h_RH.Fill(timing10);
            }
        }

        while (mcs_1.Next())
        {

            if (*detnum1 == 11)
            {
                timing5 = *time1 - timestart;
                if (timing5 > 0)
                {
                    v5.push_back(timing5);
                    PMT5.Fill(timing5);
                    vector_pair3.push_back(make_pair(timing5, 5.));
                }
            }

            else if (*detnum1 == 12)
            {
                timing6 = *time1 - timestart;
                if (timing6 > 0)
                {
                    v6.push_back(timing6);
                    PMT6.Fill(timing6);
                    vector_pair3.push_back(make_pair(timing6, 6.));
                }
            }
            else if (*detnum1 == 13)
            {
                timing7 = *time1 - timestart;
                if (timing7 > 0)
                {
                    v7.push_back(timing7);
                    PMT7.Fill(timing7);
                    vector_pair4.push_back(make_pair(timing7, 7.));
                }
            }
            else if (*detnum1 == 14)
            {
                timing8 = *time1 - timestart;
                if (timing8 > 0)
                {

                    v8.push_back(timing8);
                    PMT8.Fill(timing8);
                    vector_pair4.push_back(make_pair(timing8, 8.));
                }
            }
            else if (*detnum1 == 15)
            {
                timing9 = *time1 - timestart;
                h_GV.Fill(timing9);
            }
        }

        while (mcs_2.Next())
        {

            if (*detnum2 == 21)
            {
                timing10 = *time2;
                h_AC1.Fill(timing10);
            }
            else if (*detnum2 == 22)
            {

                timing11 = *time2;
                h_AC2.Fill(timing11);
            }
            else if (*detnum2 == 23)
            {
                timing12 = *time2;
                h_AC3.Fill(timing12);
            }
        }

        myFile->Close();

        std::vector<std::vector<double>> vec;

        vec.push_back(v1);
        vec.push_back(v2);
        vec.push_back(v3);
        vec.push_back(v4);
        vec.push_back(v5);
        vec.push_back(v6);
        vec.push_back(v7);
        vec.push_back(v8);

        /*Create a TTree with all the UCN information*/

        /*PMT12*/

        Coincidence vector12;

        ElectronicsDeadTime Corr_vec12;

        HistoGramGenerator pmt12;
        HistoGramGenerator pmt12dead;

        DeadTimeCorrection deadtime12;
        DeadTimeCorrection deadonly12;

        TTreeWriter tree12;
        PileUpCorrection PileUp12;

        Corr_vec12.VectorToFilter(vector_pair1);

        std::vector<std::pair<double, int>> FilterPair12 = Corr_vec12.ElectronicsDeadtimeCorrector();

        vector12.SetPairedVector(FilterPair12);

        std::vector<UCN> Coincidence12 = vector12.CoincidenceVectorPair();

        PileUp12.SetPileUpVector(Coincidence12);
        pmt12.SetVector(Coincidence12);
        pmt12dead.SetVector(Coincidence12);

        TH1F TDiff12 = Corr_vec12.TimeDifference();
        TDiff12.SetName("TDiff12");
        TH1F PMT12_event_length = pmt12.HistOutPutEventlength();

        /*Pileup+DeadTime*/
        TH1F PMT12 = pmt12.HistOutPut();
        TH1F PMT12_weighted = pmt12.HistOutPutWeighted();

        /*Deadtime only*/
        TH1F PMT12_dead = pmt12dead.HistOutDead();
        TH1F PMT12_weighted_dead = pmt12dead.HistOutPutDeadWeighted();

        PMT12_dead.SetName("PMT12_dead");
        PMT12_weighted_dead.SetName("PMT12_weighted_dead");

        PMT12.SetName("PMT12");
        PMT12_weighted.SetName("PMT12_weighted");
        PMT12_event_length.SetName("PMT12_event_length");

        TH1F PMT12_PileUpCorrected = PileUp12.PileUpCorrectedHistogram();
        PMT12_PileUpCorrected.SetName("PMT12_PileUpCorrected");
        deadtime12.SetHistograms(PMT12_PileUpCorrected, PMT12_weighted);

        TH1F PMT12_DeadCorrected = deadtime12.DeadTimeCorrectedHistogram();
        PMT12_DeadCorrected.SetName("PMT12_DeadCorrected");
        tree12.SetVector(Coincidence12, 1., 2.);

        deadonly12.SetHistograms(PMT12_dead, PMT12_weighted_dead);
        TH1F PMT12_DeadCorrectedOnly = deadonly12.DeadTimeCorrectedHistogram();
        PMT12_DeadCorrectedOnly.SetName("PMT12_DeadCorrectedOnly");

        TH1F PMT12_photon = pmt12.HistOutPhoton();
        PMT12_photon.SetName("PMT12_photon");

        /*PMT34*/

        Coincidence vector34;

        ElectronicsDeadTime Corr_vec34;

        HistoGramGenerator pmt34;
        HistoGramGenerator pmt34dead;

        DeadTimeCorrection deadtime34;
        DeadTimeCorrection deadonly34;

        TTreeWriter tree34;
        PileUpCorrection PileUp34;

        Corr_vec34.VectorToFilter(vector_pair2);

        std::vector<std::pair<double, int>> FilterPair34 = Corr_vec34.ElectronicsDeadtimeCorrector();

        vector34.SetPairedVector(FilterPair34);

        std::vector<UCN> Coincidence34 = vector34.CoincidenceVectorPair();

        PileUp34.SetPileUpVector(Coincidence34);

        pmt34.SetVector(Coincidence34);
        pmt34dead.SetVector(Coincidence34);

        TH1F TDiff34 = Corr_vec34.TimeDifference();
        TDiff34.SetName("TDiff34");

        TH1F PMT34_event_length = pmt34.HistOutPutEventlength();

        /*Pileup+DeadTime*/
        TH1F PMT34 = pmt34.HistOutPut();
        TH1F PMT34_weighted = pmt34.HistOutPutWeighted();

        /*Deadtime only*/
        TH1F PMT34_dead = pmt34dead.HistOutDead();
        TH1F PMT34_weighted_dead = pmt34dead.HistOutPutDeadWeighted();

        PMT34_dead.SetName("PMT34_dead");
        PMT34_weighted_dead.SetName("PMT34_weighted_dead");

        PMT34.SetName("PMT34");
        PMT34_weighted.SetName("PMT34_weighted");
        PMT34_event_length.SetName("PMT34_event_length");

        TH1F PMT34_PileUpCorrected = PileUp34.PileUpCorrectedHistogram();
        PMT34_PileUpCorrected.SetName("PMT34_PileUpCorrected");
        deadtime34.SetHistograms(PMT34_PileUpCorrected, PMT34_weighted);

        TH1F PMT34_DeadCorrected = deadtime34.DeadTimeCorrectedHistogram();
        PMT34_DeadCorrected.SetName("PMT34_DeadCorrected");
        tree34.SetVector(Coincidence34, 3., 4.);

        deadonly34.SetHistograms(PMT34_dead, PMT34_weighted_dead);
        TH1F PMT34_DeadCorrectedOnly = deadonly34.DeadTimeCorrectedHistogram();
        PMT34_DeadCorrectedOnly.SetName("PMT34_DeadCorrectedOnly");

        TH1F PMT34_photon = pmt34.HistOutPhoton();
        PMT34_photon.SetName("PMT34_photon");

        /*PMT56*/

        Coincidence vector56;
        ElectronicsDeadTime Corr_vec56;
        HistoGramGenerator pmt56;
        HistoGramGenerator pmt56dead;
        DeadTimeCorrection deadtime56;
        DeadTimeCorrection deadonly56;
        TTreeWriter tree56;
        PileUpCorrection PileUp56;

        Corr_vec56.VectorToFilter(vector_pair3);

        std::vector<std::pair<double, int>> FilterPair56 = Corr_vec56.ElectronicsDeadtimeCorrector();

        vector56.SetPairedVector(FilterPair56);

        std::vector<UCN> Coincidence56 = vector56.CoincidenceVectorPair();

        PileUp56.SetPileUpVector(Coincidence56);

        pmt56.SetVector(Coincidence56);
        pmt56dead.SetVector(Coincidence56);

        TH1F TDiff56 = Corr_vec56.TimeDifference();
        TDiff56.SetName("TDiff56");

        TH1F PMT56_event_length = pmt56.HistOutPutEventlength();

        /*Pileup+DeadTime*/
        TH1F PMT56 = pmt56.HistOutPut();
        TH1F PMT56_weighted = pmt56.HistOutPutWeighted();

        /*Deadtime only*/
        TH1F PMT56_dead = pmt56dead.HistOutDead();
        TH1F PMT56_weighted_dead = pmt56dead.HistOutPutDeadWeighted();

        PMT56_dead.SetName("PMT56_dead");
        PMT56_weighted_dead.SetName("PMT56_weighted_dead");

        PMT56.SetName("PMT56");
        PMT56_weighted.SetName("PMT56_weighted");
        PMT56_event_length.SetName("PMT56_event_length");

        TH1F PMT56_PileUpCorrected = PileUp56.PileUpCorrectedHistogram();
        PMT56_PileUpCorrected.SetName("PMT56_PileUpCorrected");
        deadtime56.SetHistograms(PMT56_PileUpCorrected, PMT56_weighted);

        TH1F PMT56_DeadCorrected = deadtime56.DeadTimeCorrectedHistogram();
        PMT56_DeadCorrected.SetName("PMT56_DeadCorrected");
        tree56.SetVector(Coincidence56, 5., 6.);

        deadonly56.SetHistograms(PMT56_dead, PMT56_weighted_dead);
        TH1F PMT56_DeadCorrectedOnly = deadonly56.DeadTimeCorrectedHistogram();
        PMT56_DeadCorrectedOnly.SetName("PMT56_DeadCorrectedOnly");

        TH1F PMT56_photon = pmt56.HistOutPhoton();
        PMT56_photon.SetName("PMT56_photon");

        /*PMT78*/

        Coincidence vector78;

        ElectronicsDeadTime Corr_vec78;

        HistoGramGenerator pmt78;
        HistoGramGenerator pmt78dead;

        DeadTimeCorrection deadtime78;
        DeadTimeCorrection deadonly78;

        TTreeWriter tree78;
        PileUpCorrection PileUp78;

        Corr_vec78.VectorToFilter(vector_pair4);

        std::vector<std::pair<double, int>> FilterPair78 = Corr_vec78.ElectronicsDeadtimeCorrector();

        vector78.SetPairedVector(FilterPair78);

        std::vector<UCN> Coincidence78 = vector78.CoincidenceVectorPair();

        PileUp78.SetPileUpVector(Coincidence78);

        pmt78.SetVector(Coincidence78);
        pmt78dead.SetVector(Coincidence78);

        TH1F TDiff78 = Corr_vec78.TimeDifference();
        TDiff78.SetName("TDiff78");

        TH1F PMT78_event_length = pmt78.HistOutPutEventlength();

        /*Pileup+DeadTime*/
        TH1F PMT78 = pmt78.HistOutPut();
        TH1F PMT78_weighted = pmt78.HistOutPutWeighted();

        /*Deadtime only*/
        TH1F PMT78_dead = pmt78dead.HistOutDead();
        TH1F PMT78_weighted_dead = pmt78dead.HistOutPutDeadWeighted();

        PMT78_dead.SetName("PMT78_dead");
        PMT78_weighted_dead.SetName("PMT78_weighted_dead");

        PMT78.SetName("PMT78");
        PMT78_weighted.SetName("PMT78_weighted");
        PMT78_event_length.SetName("PMT78_event_length");

        TH1F PMT78_PileUpCorrected = PileUp78.PileUpCorrectedHistogram();
        PMT78_PileUpCorrected.SetName("PMT78_PileUpCorrected");
        deadtime78.SetHistograms(PMT78_PileUpCorrected, PMT78_weighted);

        TH1F PMT78_DeadCorrected = deadtime78.DeadTimeCorrectedHistogram();
        PMT78_DeadCorrected.SetName("PMT78_DeadCorrected");
        tree78.SetVector(Coincidence78, 7., 8.);

        deadonly78.SetHistograms(PMT78_dead, PMT78_weighted_dead);
        TH1F PMT78_DeadCorrectedOnly = deadonly78.DeadTimeCorrectedHistogram();
        PMT78_DeadCorrectedOnly.SetName("PMT78_DeadCorrectedOnly");

        TH1F PMT78_photon = pmt78.HistOutPhoton();
        PMT78_photon.SetName("PMT78_photon");

        TFile f(file_output.c_str(), "recreate");

        f.cd();

        PMT12.Write();
        PMT12_dead.Write();
        PMT12_weighted_dead.Write();
        PMT12_weighted.Write();
        PMT12_DeadCorrected.Write();
        PMT12_event_length.Write();
        PMT12_photon.Write();
        PMT12_PileUpCorrected.Write();
        PMT12_DeadCorrectedOnly.Write();
        tree12.TTreeOutput();
        TDiff12.Write();

        PMT34.Write();
        PMT34_dead.Write();
        PMT34_weighted_dead.Write();
        PMT34_weighted.Write();
        PMT34_DeadCorrected.Write();
        PMT34_event_length.Write();
        PMT34_photon.Write();
        PMT34_PileUpCorrected.Write();
        PMT34_DeadCorrectedOnly.Write();
        tree34.TTreeOutput();
        TDiff34.Write();

        PMT56.Write();
        PMT56_dead.Write();
        PMT56_weighted_dead.Write();
        PMT56_weighted.Write();
        PMT56_DeadCorrected.Write();
        PMT56_event_length.Write();
        PMT56_photon.Write();
        PMT56_PileUpCorrected.Write();
        PMT56_DeadCorrectedOnly.Write();
        tree56.TTreeOutput();
        TDiff56.Write();

        PMT78.Write();
        PMT78_dead.Write();
        PMT78_weighted_dead.Write();
        PMT78_weighted.Write();
        PMT78_DeadCorrected.Write();
        PMT78_event_length.Write();
        PMT78_photon.Write();
        PMT78_PileUpCorrected.Write();
        PMT78_DeadCorrectedOnly.Write();
        tree78.TTreeOutput();
        TDiff78.Write();

        PMT1.Write();
        PMT2.Write();
        PMT3.Write();
        PMT4.Write();
        PMT5.Write();
        PMT6.Write();
        PMT7.Write();
        PMT8.Write();
        h_GV.Write();
        h_RH.Write();
        h_AC1.Write();
        h_AC2.Write();
        h_AC3.Write();

        f.Close();
    }
    else
    {
        myFile->Close();
    }
}
