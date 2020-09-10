/**
 *  @file   larpandoracontent/LArControlFlow/ToolMinuitFunctions.h
 *
 *  @brief  Header file for Minuit implementation.
 *
 *  $Log: $
 */

#ifndef TOOL_MINUIT_FUNCTIONS_H
#define TOOL_MINUIT_FUNCTIONS_H 1



//----------------------------------------------------------------------------------------------------------------------------------

//PHYSICAL CONSTANTS

const double K = 0.307075; // constant K in MeV cm mol^-1
const double z = 1; // charge in e
const double Z = 18; // Atomic number Z
const double A = 39.948; // Atomic mass in g mol-1
//double M = 105.7; // Mass of heavy particle in MeV
const double m_e = 0.511; // Mass of electron in MeV
const double rho = 1.396; // Density of material in g cm^-3 (here: argon density)
const double I = 0.000188; // Ionisation energy in MeV

const double C = 5.2146;
const double a = 0.19559;
const double m = 3.0;
const double X1 = 3.0;
const double X0 = 0.2000;
const double delta0 = 0.0;

//----------------------------------------------------------------------------------------------------------------------------------

void BinHitChargeVector(lar_content::TrackDirectionTool::HitChargeVector &hitChargeVector, lar_content::TrackDirectionTool::HitChargeVector &binnedHitChargeVector)
{
  //This is always commented out in TrackDirectionTool.cc  // Jan 2020: may need it!
    float binSize = (hitChargeVector.size() > 50 ? (0.5 + (hitChargeVector.size() - 50) * 2.5/300) : 0.0);

    if (binSize <= 0.0)
    {
        binnedHitChargeVector = hitChargeVector;
        return;
    }
    else if (binSize > 4.0)
        binSize = 4.0;

    float trackLength(0.f);

    for (lar_content::TrackDirectionTool::HitCharge &hitCharge : hitChargeVector)
    {
        if (hitCharge.GetLongitudinalPosition() > trackLength)
            trackLength = hitCharge.GetLongitudinalPosition();
    }

    for (float i = binSize; i <= trackLength; i += binSize)
    {
        int nHitsBin(0);
        float meanBinPosition(0.f), meanBinWidth(0.f);
        float meanBinCharge(0.f), sumSquaredSigmas(0.f);

        for (lar_content::TrackDirectionTool::HitCharge &hitCharge : hitChargeVector)
        {
            if (hitCharge.GetLongitudinalPosition() > i)
                break;

            if (hitCharge.GetLongitudinalPosition() < i && hitCharge.GetLongitudinalPosition() >= (i - binSize))
            {
                nHitsBin++;
                meanBinPosition += hitCharge.GetLongitudinalPosition();
                meanBinWidth += hitCharge.GetHitWidth();

                meanBinCharge += hitCharge.GetCharge();
                sumSquaredSigmas += (hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
            }
        }

        if (nHitsBin == 0)
            continue;

        meanBinPosition /= nHitsBin;
        meanBinWidth /= nHitsBin;
        meanBinCharge /= nHitsBin;
        sumSquaredSigmas /= (nHitsBin * nHitsBin);

        float binUncertainty(std::sqrt(sumSquaredSigmas));
        lar_content::TrackDirectionTool::HitCharge binnedHitCharge(NULL, meanBinPosition, meanBinWidth, meanBinCharge, binUncertainty);
        if (binnedHitCharge.GetCharge() > 0.1)
            binnedHitChargeVector.push_back(binnedHitCharge);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------

double DensityCorrection(double &T, double &M)
{

  //Not used in TrackDirectionTool.cc
    double p = std::sqrt((T*T) + 2*T*M);
    double gamma = std::sqrt(1 + ((p/M) * (p/M)));
    double beta = std::sqrt(1 - 1 / (gamma*gamma));
    double X = std::log10(beta*gamma);

    if (X < X0)
        return delta0;
    else if ((X > X0) && (X < X1))
        return 2 * X * std::log(10) - C + (a * (std::pow((X1 - X), m)));
    else
        return 2 * X * std::log(10) + C;
}

//----------------------------------------------------------------------------------------------------------------------------------

double BetheBloch(double &T, double &M)
{
  //used in TrackDirectionTool.cc
    double p(std::sqrt((T*T) + 2*T*M));
    double gamma(std::sqrt(1 + ((p/M) * (p/M))));
    double beta(std::sqrt(1 - 1 / (gamma*gamma)));

    double T_max(2 * m_e * (beta*gamma*beta*gamma) / (1 + 2 * gamma * m_e / M + ((m_e/M) * (m_e/M))));
    return rho * ((K * z * z * Z) / A) * (0.5*std::log(2 * m_e * T_max * (beta*gamma*beta*gamma) / (I*I) ) - (beta*beta) - (0.5*DensityCorrection(p, M))) / (beta*beta); //in MeV/cm

}

//----------------------------------------------------------------------------------------------------------------------------------

void FillLookupTable(lar_content::TrackDirectionTool::LookupTable &lookupTable, double M)
{
  //used in TrackDirectionTool.cc to set up global variables
    std::map<int, double> lookupMap;
    std::map<double, int> reverseLookupMap;

    double currentEnergy(lookupTable.GetInitialEnergy()), binWidth(lookupTable.GetBinWidth());
    int maxBin(0);

    for (double n = 0; n < 100000; ++n)
    {
        double currentdEdx = BetheBloch(currentEnergy, M);

        if ((currentdEdx * binWidth) >= currentEnergy)
        {
            double maxRange = (n * binWidth) + (currentEnergy/currentdEdx);
            lookupTable.SetMaxRange(maxRange);
            maxBin = n;

            lookupMap.insert(std::pair<int, double>(n, 0.0));
            reverseLookupMap.insert(std::pair<double, int>(0.0, n));
            break;
        }
        else
        {
            lookupMap.insert(std::pair<int, double>(n, currentEnergy));
            reverseLookupMap.insert(std::pair<double, int>(currentEnergy, n));
        }

        currentEnergy -= (currentdEdx * binWidth);
    }

    //double maxRange(lookupTable.GetMaxRange());

    //remove redundant entries to make lookup much faster
    for (std::map<int, double>::iterator it = lookupMap.begin(); it != lookupMap.end(); it++)
    {
        double n(it->first);
        double val(it->second);
        double distanceFromMaxRange((maxBin - n) * binWidth);

        if (n == 0 || n == maxBin)
            continue;
        else
        {
            double distanceMagnitude(floor(distanceFromMaxRange/2.0));
            double samplingDistance((1 + distanceMagnitude) * binWidth);

            if (!(remainder(distanceFromMaxRange, samplingDistance) == 0.0))
            {
                lookupMap.erase(n);
                reverseLookupMap.erase(val);
            }
        }
    }

    lookupTable.SetMap(lookupMap);
    lookupTable.SetReverseMap(reverseLookupMap);

    if (lookupTable.GetMaxRange() == 0.f)
        std::cout << "Warning: the lookup table max range has not been correctly set." << std::endl;
}

//----------------------------------------------------------------------------------------------------------------------------------

double GetEnergyfromLength(lar_content::TrackDirectionTool::LookupTable &lookupTable, double &trackLength)
{
  //used in TrackDirectionTool.cc
    std::map<int, double> lookupMap(lookupTable.GetMap());
    double binWidth(lookupTable.GetBinWidth());

    if (trackLength >= lookupTable.GetMaxRange())
        return 0.5; //0 energy means infinite dE/dx
    else if (trackLength <= 1.0)
      return lookupTable.GetInitialEnergy();

    int n(std::floor(trackLength/binWidth));
    std::map<int, double>::iterator nextEntryIterator(lookupMap.upper_bound(n));
    std::map<int, double>::iterator previousEntryIterator(std::prev(nextEntryIterator , 1));

    double leftLength(previousEntryIterator->first * binWidth), rightLength(nextEntryIterator->first * binWidth);
    double leftEnergy(previousEntryIterator->second), rightEnergy(nextEntryIterator->second);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(leftEnergy - rightEnergy);

    double finalEnergy(leftEnergy - (((trackLength - leftLength)/(lengthDifference)) * (energyDifference)));

    //very small energy values leadt to huge dE/dx values: truncate
    if (finalEnergy <= 2.0)
        return 2.0;
    else
        return finalEnergy;
}

//----------------------------------------------------------------------------------------------------------------------------------

double GetLengthfromEnergy(lar_content::TrackDirectionTool::LookupTable &lookupTable, double &currentEnergy)
{
  //USed in TrackDirectionTool.cc
    std::map<double, int> reverseLookupMap(lookupTable.GetReverseMap());
    double binWidth(lookupTable.GetBinWidth());

    if (currentEnergy <= 0.0)
      return lookupTable.GetMaxRange();
    else if (currentEnergy >= lookupTable.GetInitialEnergy())
        return 0.0;

    std::map<double, int>::iterator nextEntryIterator(reverseLookupMap.upper_bound(currentEnergy));
    std::map<double, int>::iterator previousEntryIterator(std::prev(nextEntryIterator , 1));

    double upperEnergy(nextEntryIterator->first), lowerEnergy(previousEntryIterator->first);
    double leftLength(nextEntryIterator->second * binWidth), rightLength(previousEntryIterator->second * binWidth);
    double lengthDifference(rightLength - leftLength);
    double energyDifference(upperEnergy - lowerEnergy);

    return leftLength + (((upperEnergy - currentEnergy)/energyDifference) * lengthDifference);
}

//----------------------------------------------------------------------------------------------------------------------------------

//void GetForwardsChiSquared(int &, double *, double &f, double *par, int )
//{

// std::cout << par[0] << f << std::endl;
  //par[0] = EndEnergy
  //par[1] = Scale
  //par[2] = Extra
  //Key Minuit equation in TrackDirectionTool.cc
  // double Ee(par[0]), L(par[1] * globalTrackLength); //energy, length

  // lar_content::TrackDirectionTool::LookupTable lookupTable = globalMuonLookupTable;
    // double M = 105.7;  //mass

    //  double Le(GetLengthfromEnergy(lookupTable, Ee));  //length
    //  // double Ls(Le - L);
    //  double Es(GetEnergyfromLength(lookupTable, Ls));    //energy

    //  double alpha((Es - Ee)/globalTotalCharge), beta(L/globalTotalHitWidth);

    //  lar_content::TrackDirectionTool::HitChargeVector binnedHitChargeVector;

    //BinHitChargeVector(*pMinuitVector, binnedHitChargeVector);

    //  double chisquared(0.0);


    //minimise this bit---------------------------------------------------------------------------------------------------------
    //  for (lar_content::TrackDirectionTool::HitCharge hitCharge : binnedHitChargeVector) //for each hit charge in the binned vector
      //  {
      //  double L_i(Ls + (par[1] * hitCharge.GetLongitudinalPosition())); //length
      //  double E_i(GetEnergyfromLength(lookupTable, L_i));       //energy

      //  double dEdx_2D(par[2] * (beta/alpha) * BetheBloch(E_i, M));   //energy per length, calculated
      //  double ChargeOverWidth(hitCharge.GetChargeOverWidth());  //energy per length, experimental

      //    chisquared += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
	//  }
    //---------------------------------------------------------------------------------------------------------------------------

    //  f = chisquared; //value to be minimised
//}

//----------------------------------------------------------------------------------------------------------------------------------

//void GetBackwardsChiSquared(int &, double *, double &f, double *par, int)
//{
//  std::cout << par[0] << f << std::endl;
  //double Ee(par[0]), L(par[1] * globalTrackLength);

  //lar_content::TrackDirectionTool::LookupTable lookupTable = globalMuonLookupTable;
  //double M = 105.7;

  ///double Le(GetLengthfromEnergy(lookupTable, Ee));
  // double Ls(Le - L);
  // double Es(GetEnergyfromLength(lookupTable, Ls));

  // double alpha((Es - Ee)/globalTotalCharge), beta(L/globalTotalHitWidth);

  //  lar_content::TrackDirectionTool::HitChargeVector binnedHitChargeVector;

    //BinHitChargeVector(*pMinuitVector, binnedHitChargeVector);

  // double chisquared(0.0);

  //  for (lar_content::TrackDirectionTool::HitCharge hitCharge : binnedHitChargeVector)
  //  {
  //     double L_i(Ls + (par[1] * (globalTrackLength - hitCharge.GetLongitudinalPosition())));
  //     double E_i(GetEnergyfromLength(lookupTable, L_i));

  //     double dEdx_2D(par[2] * (beta/alpha) * BetheBloch(E_i, M));
  //     double ChargeOverWidth(hitCharge.GetChargeOverWidth());

  //     chisquared += ( (ChargeOverWidth - dEdx_2D) * (ChargeOverWidth - dEdx_2D) )/(hitCharge.GetUncertainty() * hitCharge.GetUncertainty());
  // }

  //  f = chisquared;
//}

//----------------------------------------------------------------------------------------------------------------------------------

#endif
