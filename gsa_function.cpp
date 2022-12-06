// Markus Buchholz

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <math.h>
#include <random>
#include <numeric>

//--------Path Planner--------------------------------------------------------------

float xmin = -5.0; // 0.0;
float xmax = 5.0;  // 50.0;
float ymin = -5.0; // 0.0;
float ymax = 5.0;  // 50.0;

float obsX = 25.0;
float obsY = 25.0;
float obsR = 5.0;

float goalX = 45.0;
float goalY = 45.0;

float startX = 2.0;
float startY = 2.0;

//----------------------------------------------------------------------------------
int EVOLUTIONS = 100;
int MASSES = 20;

float ALPHA = 20.0;
float G0 = 100.0;
float EPISLION = 0.000001;

//--------------------------------------------------------------------------------

struct Pos
{

    float x;
    float y;
};
//--------------------------------------------------------------------------------

struct Velo
{

    float x;
    float y;
};
//--------------------------------------------------------------------------------

float generateNormalRandom()
{

    std::random_device engine;
    std::mt19937 gen(engine());
    std::normal_distribution<float> distrib(-1.0, 1.0);
    return distrib(gen);
}

//--------------------------------------------------------------------------------

float euclid(Pos a, Pos b)
{

    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}
//--------------------------------------------------------------------------------

float generateRandom()
{

    std::random_device engine;
    std::uniform_real_distribution<float> distrib(0, 1.0);
    return distrib(engine);
}

//--------------------------------------------------------------------------------
float valueGenerator(float low, float high)
{

    return low + generateRandom() * (high - low);
}

//--------------------------------------------------------------------------------

std::vector<Pos> initPosXY()
{

    std::vector<Pos> pos;

    for (int ii = 0; ii < MASSES; ii++)
    {

        pos.push_back({valueGenerator(xmin, xmax), valueGenerator(ymin, ymax)});
    }

    return pos;
}

//--------------------------------------------------------------------------------

std::vector<float> function(std::vector<Pos> pos)
{
    std::vector<float> funcValue;

    for (auto &ii : pos)
    {

        funcValue.push_back(ii.x * ii.y);
    }

    return funcValue;
}

//--------------------------------------------------------------------------------

float func(Pos pos)
{

    return pos.x * pos.y;
}

//--------------------------------------------------------------------------------

Pos positionUpdateCheck(Pos actualPos)
{

    Pos Pnew = actualPos;

    if (Pnew.x < xmin)
    {
        Pnew.x = xmin;
    }

    if (Pnew.x > xmax)
    {
        Pnew.x = xmax;
    }

    if (Pnew.y < ymin)
    {
        Pnew.y = ymin;
    }

    if (Pnew.y > ymax)
    {
        Pnew.y = ymax;
    }

    return Pnew;
}

//-------------------------------------------------------------------------------
bool compareMin(std::pair<Pos, float> a, std::pair<Pos, float> b)
{

    return a.second < b.second;
}

//-------------------------------------------------------------------------------
bool compareMax(std::pair<Pos, float> a, std::pair<Pos, float> b)
{

    return a.second > b.second;
}

//-------------------------------------------------------------------------------

// min
std::tuple<Pos, float> findBestPosFuncValue(std::vector<Pos> positions, std::vector<float> func)
{

    std::vector<std::pair<Pos, float>> best;

    for (int ii = 0; ii < func.size(); ii++)
    {

        best.push_back(std::pair<Pos, float>(positions[ii], func[ii]));
    }

    std::sort(best.begin(), best.end(), compareMin);

    return best[0];
}

//--------------------------------------------------------------------------------
// max
std::tuple<Pos, float> findWorstPosFuncValue(std::vector<Pos> positions, std::vector<float> func)
{

    std::vector<std::pair<Pos, float>> best;

    for (int ii = 0; ii < func.size(); ii++)
    {

        best.push_back(std::pair<Pos, float>(positions[ii], func[ii]));
    }

    std::sort(best.begin(), best.end(), compareMax);

    return best[0];
}

//--------------------------------------------------------------------------------

Velo updateVelo(Pos pos, Pos best, Velo velo, float freq)
{

    Velo veloNew;

    veloNew.x = velo.x + (pos.x - best.x) * freq;
    veloNew.y = velo.y + (pos.y - best.y) * freq;

    return veloNew;
}

//--------------------------------------------------------------------------------
// best = min, worst = max !!!! for minimization
std::vector<float> massCalculation(std::vector<float> currentValueFunction, float best, float worst)
{

    std::vector<float> masses;
    std::vector<float> masses_i;

    for (auto &ii : currentValueFunction)
    {
        float m_i = (ii - worst) / (best - worst);
        masses_i.push_back(m_i);
    }

    float sum_masses_i = std::accumulate(masses_i.begin(), masses_i.end(), 0);

    for (auto &ii : masses_i)
    {

        masses.push_back(ii / sum_masses_i);
    }

    return masses;
}

//--------------------------------------------------------------------------------
std::vector<Pos> forcesCalculation(std::vector<Pos> currentPositions, std::vector<float> masses, float Gcurrent)
{

    std::vector<Pos> forces;

    for (int ii = 0; ii < masses.size(); ii++)
    {

        Pos force_i;
        for (int jj = 0; jj < masses.size(); jj++)
        {

            if (ii != 0)
            {

                force_i.x += masses[ii] * masses[jj] * Gcurrent * (currentPositions[ii].x - currentPositions[jj].x) / (euclid(currentPositions[ii], currentPositions[jj]) + EPISLION);
                force_i.y += masses[ii] * masses[jj] * Gcurrent * (currentPositions[ii].y - currentPositions[jj].y) / (euclid(currentPositions[ii], currentPositions[jj]) + EPISLION);
            }
        }
        forces.push_back(force_i);
    }

    return forces;
}

//--------------------------------------------------------------------------------

std::vector<Pos> accCalulation(std::vector<Pos> forces)
{

    std::vector<Pos> acc;

    for (auto &ii : forces)
    {

        acc.push_back({ii.x * G0, ii.y * G0});
    }

    return acc;
}

//--------------------------------------------------------------------------------

std::vector<Pos> updataVelocities(std::vector<Pos> currentVelocities, std::vector<Pos> acc)
{

    std::vector<Pos> velo;

    for (int ii = 0; ii < currentVelocities.size(); ii++)
    {

        Pos v;

        v.x = generateRandom() * currentVelocities[ii].x + acc[ii].x;
        v.y = generateRandom() * currentVelocities[ii].y + acc[ii].y;

        velo.push_back(v);
    }

    return velo;
}

//--------------------------------------------------------------------------------
std::vector<Pos> updatePositions(std::vector<Pos> currentPositions, std::vector<Pos> currentVelocities)
{

    std::vector<Pos> pos;

    for (int ii = 0; ii < currentPositions.size(); ii++)
    {

        Pos p;

        p.x = currentPositions[ii].x + currentVelocities[ii].x;
        p.y = currentPositions[ii].y + currentVelocities[ii].y;

        pos.push_back(positionUpdateCheck(p));
    }

    return pos;
}

//--------------------------------------------------------------------------------

void runGSA()
{

    std::vector<Pos> currentPositions = initPosXY();
    std::vector<Pos> currentVelocities(MASSES, {0.0, 0.0});
    std::vector<float> currentValueFunction(MASSES, 0);

    for (int ii = 0; ii < EVOLUTIONS; ii++)
    {
        currentValueFunction = function(currentPositions);
        std::tuple<Pos, float> bestValueFun = findBestPosFuncValue(currentPositions, currentValueFunction);
        float bestValue = std::get<1>(bestValueFun);
        std::tuple<Pos, float> worstValueFun = findWorstPosFuncValue(currentPositions, currentValueFunction);
        float worstValue = std::get<1>(worstValueFun);

        float Gcurrent = G0 * std::exp(-ALPHA * ii / EVOLUTIONS);
        std::vector<float> masses = massCalculation(currentValueFunction, bestValue, worstValue);
        std::vector<Pos> forces = forcesCalculation(currentPositions, masses, Gcurrent);
        std::vector<Pos> acc = accCalulation(forces);
        currentVelocities = updataVelocities(currentVelocities, acc);
        currentPositions = updatePositions(currentPositions, currentVelocities);
    }

    for (auto &ii : currentValueFunction)
    {
        std::cout << " f = " << ii << std::endl;
    }
}

//--------------------------------------------------------------------------------

int main()
{

    runGSA();
}