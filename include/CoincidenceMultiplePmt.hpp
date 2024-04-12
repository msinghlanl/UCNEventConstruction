#ifndef COINCIDENCEMULTIPLEPMT_HPP
#define COINCIDENCEMULTIPLEPMT_HPP
#include "TH1F.h"

struct UCN
{
    double StartTime;
    double ClusterSize;
    double EventLength;
    int Na;
    int Nb;
    double Npile;
};

class Coincidence
{
    std::vector<double> v1, v2;
    std::vector<std::pair<double, int>> vector_pair;

public:
    void SetVectorPair(std::vector<double> PMTa, std::vector<double> PMTb)
    {
        v1 = PMTa;
        v2 = PMTb;
    }

    void SetPairedVector(std::vector<std::pair<double, int>> vectorpair)
    {
        vector_pair = vectorpair;
    }

    /*Methods for the UCN*/
    std::vector<UCN> SetCoincidence();
    std::vector<UCN> CoincidenceVectorPair();

    void PrintElements();

    void FillTree();

    /*Resulant Coincidence Vector*/
    std::vector<UCN> CoincidenceVector;
    std::vector<UCN> CoincidenceVector1;

    TTree *UCNTree = new TTree("UCNTree", "UCN event tree");
};

class HistoGramGenerator
{
    std::vector<UCN> vector;

public:
    void
    SetVector(std::vector<UCN> coincidence_vector)
    {
        vector = coincidence_vector;
    }
    TH1F HistOutPut();
    TH1F HistOutPutWeighted();
    TH1F HistOutPutEventlength();
    TH1F HistOutPhoton();
    TH1F HistOutPutDeadWeighted();
    TH1F HistOutDead();
};

class TTreeWriter
{
    std::vector<UCN> vector;
    int det1, det2;

public:
    void SetVector(std::vector<UCN> coincidence_vector, int pmt1, int pmt2)
    {
        vector = coincidence_vector;
        det1 = pmt1;
        det2 = pmt2;
    }

    void TTreeOutput();
};

class DeadTimeCorrection
{
    TH1F h_finer;
    TH1F h_weighted;

public:
    void SetHistograms(TH1F H1, TH1F H2)
    {
        h_finer = H1;
        h_weighted = H2;
    }
    TH1F DeadTimeCorrectedHistogram();
};

class PileUpCorrection
{
    std::vector<UCN> vector;

public:
    void SetPileUpVector(std::vector<UCN> coincidence_vector)
    {

        vector = coincidence_vector;
    }
    std::vector<double> PileUpCorrectedVector();
    TH1F PileUpCorrectedHistogram();
};

class ElectronicsDeadTime
{
    std::vector<std::pair<double, int>> pmt_singles;

public:
    void VectorToFilter(std::vector<std::pair<double, int>> input_vector)
    {
        pmt_singles = input_vector;
    }
    std::vector<std::pair<double, int>> ElectronicsDeadtimeCorrector();

    TH1F TimeDifference();
};

#endif

/*for (int bin = 1; bin <= PMT12_finer->GetNbinsX(); bin++)
   {
       double counts = PMT12_finer->GetBinContent(bin);
       double weighted_counts = h_weighted12->GetBinContent(bin);
       double length = h_event12->GetBinContent(bin);
       double correction_factor = 0.1 / (0.1 - weighted_counts);
       double final_factor = 0.;
       if (correction_factor < 1.0 || correction_factor > 1.1)
       {
           final_factor = 1.0;
       }

       else
       {
           final_factor = correction_factor;
       }
       double corrected_counts = counts * final_factor;
       PMT12_deadtime_corrected->SetBinContent(bin, corrected_counts);
   }*/