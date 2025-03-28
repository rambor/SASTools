//   Copyright 2020 Robert P. Rambo
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include <cfloat>
#include <random>
#include "utils.h"
#include "vector3.h"

using namespace std;

static const boost::regex columnAtomBranch("([1-9PABGDEZH]|\\s)");
static const boost::regex columnNumBranch("([1-9])");
static const boost::regex columnOxygenBranch("([ABGDENXH])");
static const boost::regex isWater("HOH");

// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;


/*
 * reset qvalues vector using a specified percentage of original data
 */
void sub_select(std::vector <float> & qvalues, std::vector<std::vector <float> > & iobs, float subset_fraction, int * qvaluesSize, int lmax, float qmin, float qmax){

    int * tempselected = NULL;
    int totalqValues = *qvaluesSize;
    int * tempselectedIndex = NULL;
    float * tempq = NULL;
    float incrq = 0.0, endq, startq;
    int binSize, select = 0, perBin, selectLocation = 0, newSize, atq = 0;
    int count = 0;
    int selectedIndex;

    srand( unsigned ( time(NULL) ) );

    std::vector<float> oldqvalues(totalqValues);

    for(int i=0; i<totalqValues; i++){
        oldqvalues[i] = qvalues[i];
    }

    binSize = (int)(totalqValues/lmax);
    startq = qmin;
    incrq = (qmax - qmin)/lmax;

    perBin = (int)(totalqValues/lmax*subset_fraction);
    newSize = (int)(totalqValues*subset_fraction);

    tempq = new float[newSize];
    tempselectedIndex = new int[newSize];

    if (perBin < 1) {
        perBin = 1;
    }

    for (int l = 0; l < lmax; l++){
        //
        // grab all values in current bin given by l
        //
        tempselected = new int[binSize];
        endq = startq + incrq; // *(l+1);

        // grab all the values <= startq and < endq
        // start loop from last true statement (qvalues[q] < endq)
        for (int q=atq; q < totalqValues; q++){

            if ( (qvalues[q] >= startq) && (qvalues[q] < endq) && (select < binSize) ){
                tempselected[select] = q; //store integer location of qvalue of interest
                select++;
            }
        }

        atq=atq+select;

        // randomly select from tempselected and put into selectedQ
        random_shuffle( &tempselected[0], &tempselected[select - 1], p_myrandom );

        select = 0;
        // grab required elements per bin
        for (int b=0; b < perBin; b++){
            // some indices may not be appropriate due to the binSize
            // what if binSize is 70 but only 60 datapoints meet criteria?
            if ((tempselected[b] > 0) && (tempselected[b] < atq)){
                tempselectedIndex[selectLocation] = tempselected[b];
                selectLocation++;
            }
        }

        delete[] tempselected;
        tempselected = NULL;
        startq = endq;
    }

    std::sort(&tempselectedIndex[0], &tempselectedIndex[0] + newSize);    // holds indices

    // copy values to temp2d

    vector < vector <float> > temp2d(newSize, vector<float>(2));

    for (int b=0; b < newSize; b++){

        selectedIndex = tempselectedIndex[b];

        if ((selectedIndex < totalqValues) && (selectedIndex >= 0) ) {
            qvalues[count] = oldqvalues[selectedIndex];
            temp2d[count] = iobs[selectedIndex];
            count++;
        }
    }

    qvalues.resize(count);
    temp2d.resize(count);
    delete[] tempq;
    delete[] tempselectedIndex;
    tempq = NULL;
    tempselectedIndex = NULL;

    *qvaluesSize = count;
    iobs = temp2d;

}

/*
 * Convert cartesian to spherical coordinates
 * Return 3 element array
 */
float * xyz_to_rtp(float const &tempx, float const &tempy, float const &tempz)
{
    float x, y, z, r, theta, phi;
    static float tempArray[3];
    x = tempx;
    y = tempy;
    z = tempz;

    r = sqrt(x*x + y*y + z*z);

    theta = acosf(z/r);
    if (x == 0){
        phi = 0.0;
    } else {
        phi = atan2(y, x);
    }
    tempArray[0] = r;
    tempArray[1] = theta;
    tempArray[2] = phi;

    return tempArray;
}



/*
 * Center coordinates and determine dmax
 * Calculate dmax from PDB, pass in point to the Array containing vectors of x,y,z, number of atoms
 * and the pointer to hold the center of the coordinates
 */
void dmaxFromPDB(std::vector <float> & x,
                 std::vector <float> & y,
                 std::vector <float> & z,
                 float * dmax,
                 float * centeredX, // export centered coordinates
                 float * centeredY,
                 float * centeredZ,
                 vector3 * centeringVec,
                 size_t const numAtoms
) {


    // sum each array and divide by size to get center of coordinates
    float sumX = std::accumulate(x.begin(), x.end(), 0.0f);
    float sumY = std::accumulate(y.begin(), y.end(), 0.0f);
    float sumZ = std::accumulate(z.begin(), z.end(), 0.0f);

    float centeringX = sumX/(float)numAtoms;
    float centeringY = sumY/(float)numAtoms;
    float centeringZ = sumZ/(float)numAtoms;
    *centeringVec = vector3(centeringX, centeringY, centeringZ);

    std::vector<vector3> minmax = std::vector<vector3>(14);
    std::vector<float> abcdefgh = std::vector<float>(8);
    std::fill(abcdefgh.begin(), abcdefgh.end(), 0.0);

    float * pX = x.data();
    float * pY = y.data();
    float * pZ = z.data();
    float pXValue, pYValue, pZValue;


    for (unsigned int n=0; n < numAtoms; n++) {
        //populate array
        centeredX[n] = pX[n] - centeringX;
        centeredY[n] = pY[n] - centeringY;
        centeredZ[n] = pZ[n] - centeringZ;

        pXValue = centeredX[n];
        pYValue = centeredY[n];
        pZValue = centeredZ[n];

        // test if extreme points for dmax calculation
        if (pXValue > minmax[0].x){
            minmax[0].x = pXValue;
            minmax[0].y = pYValue;
            minmax[0].z = pZValue;
        }

        if (pXValue < minmax[1].x){
            minmax[1].x = pXValue;
            minmax[1].y = pYValue;
            minmax[1].z = pZValue;
        }

        if (pYValue > minmax[2].y){
            minmax[2].x = pXValue;
            minmax[2].y = pYValue;
            minmax[2].z = pZValue;
        }

        if (pYValue < minmax[3].y){
            minmax[3].x = pXValue;
            minmax[3].y = pYValue;
            minmax[3].z = pZValue;
        }

        if (pZValue > minmax[4].z){
            minmax[4].x = pXValue;
            minmax[4].y = pYValue;
            minmax[4].z = pZValue;
        }

        if (pZValue > minmax[5].z){
            minmax[5].x = pXValue;
            minmax[5].y = pYValue;
            minmax[5].z = pZValue;
        }


        if ((pXValue - pYValue + pZValue) > abcdefgh[0]){
            minmax[6].x = pXValue;
            minmax[6].y = pYValue;
            minmax[6].z = pZValue;
            abcdefgh[0] = (pXValue - pYValue + pZValue);
        }

        if ((pXValue - pYValue - pZValue) > abcdefgh[1]){
            minmax[7].x = pXValue;
            minmax[7].y = pYValue;
            minmax[7].z = pZValue;
            abcdefgh[1] = (pXValue - pYValue - pZValue);
        }


        if ((pXValue + pYValue + pZValue) > abcdefgh[2]){
            minmax[8].x = pXValue;
            minmax[8].y = pYValue;
            minmax[8].z = pZValue;
            abcdefgh[2] = (pXValue + pYValue + pZValue);
        }

        if ((pXValue + pYValue - pZValue) > abcdefgh[3]){
            minmax[9].x = pXValue;
            minmax[9].y = pYValue;
            minmax[9].z = pZValue;
            abcdefgh[3] = (pXValue + pYValue - pZValue);
        }


        // minimizes
        if ((pXValue - pYValue + pZValue) < abcdefgh[4]){
            minmax[10].x = pXValue;
            minmax[10].y = pYValue;
            minmax[10].z = pZValue;
            abcdefgh[4] = (pXValue - pYValue + pZValue);
        }

        if ((pXValue - pYValue - pZValue) < abcdefgh[5]){
            minmax[11].x = pXValue;
            minmax[11].y = pYValue;
            minmax[11].z = pZValue;
            abcdefgh[5] = (pXValue - pYValue - pZValue);
        }

        if ((pXValue + pYValue + pZValue) < abcdefgh[6]) {
            minmax[12].x = pXValue;
            minmax[12].y = pYValue;
            minmax[12].z = pZValue;
            abcdefgh[6] = (pXValue + pYValue + pZValue);
        }

        if ((pXValue + pYValue - pZValue) < abcdefgh[7]) {
            minmax[13].x = pXValue;
            minmax[13].y = pYValue;
            minmax[13].z = pZValue;
            abcdefgh[7] = (pXValue + pYValue - pZValue);
        }
    }

    // calculate dmax
    *dmax = FLT_MIN;
    float tdmax;
    for(unsigned int i =0 ; i<14; i++){
        vector3 * pVec1 = &minmax[i];
        for(unsigned int j=(i+1) ; j<14; j++){
            tdmax =(*pVec1 - minmax[j]).length();
            if (tdmax > *dmax){
                *dmax = tdmax;
            }
        }
    }

    pX = nullptr;
    pY = nullptr;
    pZ = nullptr;
//    std::cout << "\td_max estimated as: " << *dmax << std::endl;
//    std::cout << "\tCentering PDB model at (" << centerX << ", " << centerY  << ", " << centerZ << ")"<< std::endl;
}



/*
 * Calculate water lattice, return 5 element vector containing waters (x,y,z, nieghboring atom and residue type)
 *
 */

void waterLattice(
        float const * const x,
        float const * const y,
        float const * const z,
        float lowerBound,
        float upperBound,
        float waters[][5],
//                  vector <vector<float> > &waters,
        size_t const numAtoms,
        int * numWaters)
{
    std::vector<double>::iterator it;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> randomDis(0.0, 1);

//	float fractH, fractK, fractL, fractHtest, fractKtest, fractLtest, h, k, l, testH, testK, testL, ht, kt, lt;
    float h=0, k=0, l=0, ht=0, kt=0, lt=0;

    float minimum;
    float maxH, minH, maxK, minK, maxL, minL, garbage;

    minH = *min_element(x, x+numAtoms);
    maxH = *max_element(x, x+numAtoms);
    minK = *min_element(y, y+numAtoms);
    maxK = *max_element(y, y+numAtoms);
    minL = *min_element(z, z+numAtoms);
    maxL = *max_element(z, z+numAtoms);

    const float boxEdge = 0.0694444f; //  = 1/14.4

    garbage = modf((*min_element(x, x+numAtoms)-14.4f)*boxEdge, &minH); // determine hkl limits, i.e., max and min
    garbage = modf((*max_element(x, x+numAtoms)+14.4f)*boxEdge, &maxH);
    garbage = modf((*min_element(y, y+numAtoms)-14.4f)*boxEdge, &minK);
    garbage = modf((*max_element(y, y+numAtoms)+14.4f)*boxEdge, &maxK);
    garbage = modf((*min_element(z, z+numAtoms)-14.4f)*boxEdge, &minL);
    garbage = modf((*max_element(z, z+numAtoms)+14.4f)*boxEdge, &maxL);


    std::list<int>::const_iterator iter;
    //declare the starting grid of TIP3 waters at 25 degrees in the 14.4 A box centered at (0,0,0)
    double wxArr[] = {-4.313255814,3.610744186,-5.402255814,3.362744186,3.862744186,1.828744186,2.871744186,6.827744186,1.759744186,2.014744186,3.218744186,-2.445255814,-3.660255814,-3.944255814,-1.605255814,-4.719255814,-4.353255814,-6.316255814,-0.036255814,-3.035255814,5.429744186,4.113744186,6.070744186,0.865744186,-3.743255814,-2.524255814,5.729744186,-0.396255814,-2.409255814,-1.575255814,-3.377255814,2.966744186,3.486744186,6.004744186,-1.613255814,-0.206255814,5.049744186,-4.642255814,-1.922255814,-4.554255814,-1.602255814,4.716744186,6.346744186,0.250744186,-6.774255814,5.376744186,-0.095255814,6.469744186,3.899744186,0.465744186,4.136744186,-2.179255814,6.803744186,-3.889255814,5.209744186,-3.946255814,2.986744186,0.817744186,1.495744186,5.975744186,1.336744186,1.952744186,-0.624255814,4.252744186,3.343744186,6.972744186,1.836744186,5.350744186,-3.983255814,-0.808255814,-0.923255814,-1.913255814,-4.884255814,-6.375255814,-6.614255814,-6.150255814,-4.785255814,-3.625255814,-5.501255814,-5.749255814,0.508744186,-5.930255814,-1.784255814,-6.388255814,1.514744186,0.223744186};
    double wyArr[] = {1.637488372,1.702488372,-0.525511628,0.397488372,-1.841511628,5.770488372,5.669488372,-5.287511628,3.374488372,4.760488372,-6.566511628,0.918488372,-2.519511628,3.118488372,5.890488372,3.364488372,-1.021511628,-2.737511628,-4.982511628,0.492488372,2.244488372,2.403488372,-2.671511628,2.394488372,-5.767511628,-4.601511628,6.157488372,0.340488372,-1.489511628,7.048488372,1.146488372,4.234488372,2.337488372,-4.292511628,1.860488372,3.931488372,-6.285511628,-3.431511628,-5.595511628,-3.701511628,-0.452511628,6.302488372,6.589488372,-6.978511628,3.331488372,5.200488372,4.422488372,3.343488372,-5.469511628,-2.698511628,-3.055511628,-3.125511628,-2.582511628,3.952488372,-1.782511628,4.612488372,-1.601511628,-3.688511628,-3.691511628,3.304488372,-1.654511628,-3.120511628,-2.062511628,-4.933511628,-0.557511628,0.523488372,0.103488372,1.797488372,-6.257511628,2.595488372,-3.537511628,4.478488372,6.257488372,-4.425511628,0.159488372,5.113488372,-3.284511628,0.892488372,-0.670511628,2.178488372,-6.786511628,-0.126511628,-1.451511628,4.367488372,-6.374511628,2.965488372};
    double wzArr[] = {2.381918605,-2.261081395,-5.304081395,-4.843081395,-6.229081395,-1.926081395,2.943918605,-4.255081395,1.579918605,5.486918605,-5.194081395,0.529918605,2.567918605,5.910918605,-2.658081395,-4.623081395,0.292918605,3.331918605,6.163918605,5.528918605,-5.032081395,0.680918605,-0.688081395,-2.925081395,-3.554081395,1.047918605,4.920918605,2.398918605,-2.224081395,5.964918605,-3.846081395,-3.813081395,4.067918605,2.901918605,-1.910081395,-0.418081395,-2.305081395,-2.189081395,-6.748081395,5.662918605,-5.304081395,-6.264081395,-4.041081395,1.904918605,6.220918605,-1.877081395,-4.894081395,-0.087081395,4.542918605,-5.851081395,-3.776081395,6.465918605,-3.293081395,-1.293081395,3.156918605,2.822918605,-1.489081395,-1.547081395,4.287918605,4.529918605,6.309918605,1.566918605,0.207918605,0.708918605,1.640918605,-0.603081395,4.042918605,6.661918605,4.540918605,-6.883081395,-3.543081395,4.746918605,-0.070081395,-0.133081395,1.892918605,-6.310081395,-6.047081395,-6.931081395,4.991918605,-0.074081395,-3.769081395,-2.169081395,4.209918605,3.704918605,-0.896081395,6.567918605};
    int numberOfWaters = sizeof(wxArr)/sizeof(double);
    int nPosition = 0;

    std::vector<float> wVx(wxArr, wxArr+numberOfWaters); //make vector from the water coordinate array
    std::vector<float> wVy(wyArr, wyArr+numberOfWaters);
    std::vector<float> wVz(wzArr, wzArr+numberOfWaters);

    transform (wVx.begin(), wVx.end(), wVx.begin(), std::bind2nd( std::plus<double>(), 7.2)); //shift box such that the corner is (0,0,0)
    transform (wVy.begin(), wVy.end(), wVy.begin(), std::bind2nd( std::plus<double>(), 7.2));
    transform (wVz.begin(), wVz.end(), wVz.begin(), std::bind2nd( std::plus<double>(), 7.2));

    std::vector<float> tempVh_x(0); // temp - holds the atom coordinates that correspond to the working (h, k, l)
    std::vector<float> tempVk_y(0);
    std::vector<float> tempVl_z(0);
    std::vector<int> usedAtoms(0);
    std::vector<double> distances(0);

    ht = h*14.4f;
    kt = k*14.4f;
    lt = l*14.4f;

    // initialize random seed
    float test = 0.0f;


    for (int ch = (int) minH; ch <= (int) maxH; ch++) {

        for (int ck = (int) minK; ck <= (int) maxK; ck++) {

            for (int cl = (int) minL; cl <= (int) maxL; cl++) {

                float wx, wy, wz, tx, ty, tz;
                // go through each water and determine if it is acceptable


                for(int n = 0; n < numberOfWaters; n++) {

                    distances.clear();
                    wx = wVx[n] + ch*14.4f;
                    wy = wVy[n] + ck*14.4f;
                    wz = wVz[n] + cl*14.4f;


                    for (int atom = 0; atom < numAtoms; atom++) {
                        // calculate the distances from select water to all atoms
                        tx = x[atom];
                        ty = y[atom];
                        tz = z[atom];

                        distances.push_back(  (tx - wx)*(tx - wx)
                                              + (ty - wy)*(ty - wy)
                                              + (tz - wz)*(tz - wz) );

                    }

                    // determine if minimum distance is greater than lower and less than upper bounds
                    minimum = *(min_element(distances.begin(), distances.end()));


                    if (minimum > lowerBound && minimum < upperBound) {
                        it = find(distances.begin(), distances.end(), minimum);
                        // Need to identify the protein/RNA atom that is closest to the accepted water
                        nPosition = distance(distances.begin(), it);

                        waters[*numWaters][0] = wx;
                        waters[*numWaters][1] = wy;
                        waters[*numWaters][2] = wz;
                        //
                        // lowering the cut-off pushes more waters into diffuse bulk
                        if (minimum > 7.25) { // 2.8^2 = 7.84; 2.69^2 = 7.25, 2.6 = 6.76, 2.5^2 = 6.25, 2.4 = 5.76
                            waters[*numWaters][3] = -1.0f;
                        } else {
                            waters[*numWaters][3] = nPosition;
                        }

                        *numWaters = *numWaters + 1;

                    } else if ((minimum > upperBound) && (minimum < 12.25)) {
                        // 12.25 = 3.5*3.5
                        // accept with some probability proporitional to distance
                        test =  expf(-(minimum - upperBound)/(upperBound));

                        if (randomDis(gen) < test) {
//                          cout << "Prob " <<  acceptance << " " << test << " " << minimum << " " << upperBound <<  endl;
//
//                            waters[*numWaters][0] = wx;
//                            waters[*numWaters][1] = wy;
//                            waters[*numWaters][2] = wz;
//                            waters[*numWaters][3] = -1.0;
//
//                            *numWaters = *numWaters + 1;
//
                        }
                    }
                }
            } // end L
        }  // end K
    } // end H
}

float assignOccupancy ( std::string * neighboringAtom, std::string * neighboringResi){
    float occupancy = 0;
    string atom, residue;
//	boost::algorithm::trim(*neighboringAtom);
//	boost::algorithm::trim(*neighboringResi);

    //RNA bases
    if ((*neighboringResi == "rA")||(*neighboringResi == "rG") || (*neighboringResi == "rC") || (*neighboringResi == "rU") || \
        (*neighboringResi == "ADE") || (*neighboringResi == "GUA") || (*neighboringResi == "CYT") || (*neighboringResi == "URI") ) {
        // Adapted from Hydration Sites of unpair RNA bases (BMC Structural Biology 2011 11:41)
        if (*neighboringAtom == "N1") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C2") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "O2") {
            occupancy = 0.2;
        } else if (*neighboringAtom == "N2") {
            occupancy = 0.185;
        } else if (*neighboringAtom == "N3") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C4") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "O4") {
            occupancy = 0.36;
        } else if (*neighboringAtom == "C5") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C6") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "N6") {
            occupancy = 0.65;
        } else if (*neighboringAtom == "N7") {
            occupancy = 0.65;
        } else if (*neighboringAtom == "C8") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "N9") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C1\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C2\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C3\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C4\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "O2\'") {
            occupancy = 0.348;
        } else if (*neighboringAtom == "O3\'") {
            occupancy = 0.5;
        } else if (*neighboringAtom == "O4\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "C5\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "O5\'") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "P") {
            occupancy = 0.0;
        } else if (*neighboringAtom == "O1P") {
            occupancy = 0.5;
        } else if (*neighboringAtom == "O2P") {
            occupancy = 0.5;
        } else {
            occupancy = 1.0;
        }
    } else if ((*neighboringResi == "dA")||(*neighboringResi == "dG") || (*neighboringResi == "dC") || (*neighboringResi == "dT")) {

        if (*neighboringAtom == "N1") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C2") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O2") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "N2") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "N3") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C4") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O4") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C5") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C6") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "N6") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "N7") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C8") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "N9") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C1\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C2\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C3\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C4\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O2\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O3\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O4\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "C5\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O5\'") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "P") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O1P") {
            occupancy = 1.0;
        } else if (*neighboringAtom == "O2P") {
            occupancy = 1.0;
        } else {
            occupancy = 1.0;
        }

    } else if (*neighboringResi == "ALA") {
        if (*neighboringAtom == "N") {
            occupancy = 0.355;
//			occupancy = 0.595819;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.072;
//			occupancy = 0.268328;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.007;
//			occupancy = 0.083666;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.584;
//			occupancy = 0.764198;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.232;
//			occupancy = 0.4816637;
        }
    } else if (*neighboringResi == "ARG") {
        if (*neighboringAtom == "N") {
            occupancy = 0.470;
//			occupancy = 0.685565;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.039;
//			occupancy = 0.197484;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.012;
//			occupancy = 0.109544;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.606;
//			occupancy = 0.77846;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.033;
//			occupancy = 0.181659;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.085;
//			occupancy = 0.291547;
        } else if (*neighboringAtom == "CD") {
            occupancy = 0.172;
//			occupancy = 0.4147288;
        } else if (*neighboringAtom == "NE") {
            occupancy = 0.237;
//			occupancy = 0.486826;
        } else if (*neighboringAtom == "CZ") {
            occupancy = 0.026;
//			occupancy = 0.161245;
        } else if (*neighboringAtom == "NH1") {
            occupancy = 0.440;
//			occupancy = 0.663324;
        } else if (*neighboringAtom == "NH2") {
            occupancy = 0.445;
//			occupancy = 0.667;
        }
    } else if (*neighboringResi == "ASN") {
        if (*neighboringAtom == "N") {
            occupancy = 0.367;
//            occupancy = 0.605805;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.021;
//            occupancy = 0.1449137;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.017;
//            occupancy = 0.130384;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.448;
//            occupancy = 0.669328;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.088;
//            occupancy = 0.296647;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.040;
//            occupancy = 0.2;
        } else if (*neighboringAtom == "OD1") {
            occupancy = 0.459;
//            occupancy = 0.677495;
        } else if (*neighboringAtom == "ND2") {
            occupancy = 0.464;
//            occupancy = 0.6811754;
        }
    } else if (*neighboringResi == "ASP") {
        if (*neighboringAtom == "N") {
            occupancy = 0.301;
//            occupancy = 0.5488346;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.011;
//            occupancy = 0.10488;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.006;
//            occupancy = 0.077459;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.503;
//            occupancy = 0.709225;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.082;
//            occupancy = 0.286356;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.037;
//            occupancy = 0.192353;
        } else if (*neighboringAtom == "OD1") {
            occupancy = 0.546;
//            occupancy = 0.738918;
        } else if (*neighboringAtom == "OD2") {
            occupancy = 0.546;
//            occupancy = 0.738918;
        }
    } else if (*neighboringResi == "CYS" || *neighboringResi == "CYX") {
        if (*neighboringAtom == "N") {
            occupancy = 0.400;
//            occupancy = 0.632455;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.036;
//            occupancy = 0.189736;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.014;
//            occupancy = 0.118321;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.583;
//            occupancy = 0.763544;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.100;
//            occupancy = 0.31622;
        } else if (*neighboringAtom == "SG") {
            occupancy = 0.139;
//            occupancy = 0.372827;
        }
    } else if (*neighboringResi == "GLN") {
        if (*neighboringAtom == "N") {
            occupancy = 0.415;
//            occupancy = 0.644205;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.040;
//            occupancy = 0.2;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.009;
//            occupancy = 0.094868;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.502;
//            occupancy = 0.708519;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.064;
//            occupancy = 0.25298;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.116;
//            occupancy = 0.340589;
        } else if (*neighboringAtom == "CD") {
            occupancy = 0.020;
//            occupancy = 0.141421;
        } else if (*neighboringAtom == "OE1") {
            occupancy = 0.453;
//            occupancy = 0.67305;
        } else if (*neighboringAtom == "NE2") {
            occupancy = 0.416;
//            occupancy = 0.6449806;
        }
    } else if (*neighboringResi == "GLU") {
        if (*neighboringAtom == "N") {
            occupancy = 0.343;
//            occupancy = 0.585662;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.040;
//            occupancy = 0.2;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.006;
//            occupancy = 0.077459;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.513;
//            occupancy = 0.7162;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.070;
//            occupancy = 0.264575;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.109;
//            occupancy = 0.330151;
        } else if (*neighboringAtom == "CD") {
            occupancy = 0.014;
//            occupancy = 0.11832;
        } else if (*neighboringAtom == "OE1") {
            occupancy = 0.481;
//            occupancy = 0.693541;
        } else if (*neighboringAtom == "OE2") {
            occupancy = 0.469;
//            occupancy = 0.684836;
        }
    } else if (*neighboringResi == "GLY") {
        if (*neighboringAtom == "N") {
            occupancy = 0.312;
//            occupancy = 0.558569;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.139;
//            occupancy = 0.372827;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.028;
//            occupancy = 0.167332;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.521;
//            occupancy = 0.721803;
        }
    } else if (*neighboringResi == "HIS" || *neighboringResi == "HID" || *neighboringResi == "HID" || *neighboringResi == "HIP") {
        if (*neighboringAtom == "N") {
            occupancy = 0.359;
//            occupancy = 0.599166;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.023;
//            occupancy = 0.15165;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.022;
//            occupancy = 0.148323;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.442;
//            occupancy = 0.664830;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.099;
//            occupancy = 0.314642;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.060;
//            occupancy = 0.244948;
        } else if (*neighboringAtom == "CD2") {
            occupancy = 0.140;
//            occupancy = 0.374165;
        } else if (*neighboringAtom == "ND1") {
            occupancy = 0.441;
//            occupancy = 0.664078;
        } else if (*neighboringAtom == "CE1") {
            occupancy = 0.165;
//            occupancy = 0.406202;
        } else if (*neighboringAtom == "NE2") {
            occupancy = 0.353;
//            occupancy = 0.594138;
        }
    } else if (*neighboringResi == "ILE") {
        if (*neighboringAtom == "N") {
            occupancy = 0.380;
//            occupancy = 0.61644;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.024;
//            occupancy = 0.154919;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.545;
//            occupancy = 0.73824;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.018;
//            occupancy = 0.13416;
        } else if (*neighboringAtom == "CG1") {
            occupancy = 0.031;
//            occupancy = 0.176068;
        } else if (*neighboringAtom == "CG2") {
            occupancy = 0.124;
//            occupancy = 0.352136;
        } else if (*neighboringAtom == "CD1") {
            occupancy = 0.087;
//            occupancy = 0.294957;
        }
    } else if (*neighboringResi == "LEU") {
        if (*neighboringAtom == "N") {
            occupancy = 0.384;
//            occupancy = 0.619677;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.010;
//            occupancy = 0.1;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.548;
//            occupancy = 0.74027;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.035;
//            occupancy = 0.187082;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.020;
//            occupancy = 0.141421;
        } else if (*neighboringAtom == "CD1") {
            occupancy = 0.125;
//            occupancy = 0.35355;
        } else if (*neighboringAtom == "CD2") {
            occupancy = 0.117;
//            occupancy = 0.342052;
        }
    } else if (*neighboringResi == "LYS") {
        if (*neighboringAtom == "N") {
            occupancy = 0.271;
//            occupancy = 0.520576;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.012;
//            occupancy = 0.109544;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.009;
//            occupancy = 0.094868;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.474;
//            occupancy = 0.688476;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.060;
//            occupancy = 0.2449489;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.076;
//            occupancy = 0.27568;
        } else if (*neighboringAtom == "CD") {
            occupancy = 0.081;
//            occupancy = 0.284604;
        } else if (*neighboringAtom == "CE") {
            occupancy = 0.129;
//            occupancy = 0.359165;
        } else if (*neighboringAtom == "NZ") {
            occupancy = 0.476;
//            occupancy = 0.689927;
        }
    } else if (*neighboringResi == "MSE") {
        if (*neighboringAtom == "N") {
            occupancy = 0.486;
//            occupancy = 0.697137;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.030;
//            occupancy = 0.173205;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.030;
//            occupancy = 0.173205;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.591;
//            occupancy = 0.768765;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.132;
//            occupancy = 0.3633180;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.095;
//            occupancy = 0.308220;
        } else if (*neighboringAtom == "SE") {
            occupancy = 0.015;
//            occupancy = 0.122474;
        } else if (*neighboringAtom == "CE") {
            occupancy = 0.208;
//            occupancy = 0.45607017;
        }
    } else if (*neighboringResi == "MET") {
        if (*neighboringAtom == "N") {
            occupancy = 0.486;
//            occupancy = 0.697137;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.030;
//            occupancy = 0.173205;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.030;
//            occupancy = 0.173205;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.591;
//            occupancy = 0.768765;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.132;
//            occupancy = 0.3633180;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.095;
//            occupancy = 0.308220;
        } else if (*neighboringAtom == "SD") {
            occupancy = 0.015;
//            occupancy = 0.122474;
        } else if (*neighboringAtom == "CE") {
            occupancy = 0.208;
//            occupancy = 0.45607017;
        }
    } else if (*neighboringResi == "PHE") {
        if (*neighboringAtom == "N") {
            occupancy = 0.406;
//            occupancy = 0.637181;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.014;
//            occupancy = 0.11832;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.503;
//            occupancy = 0.709224;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.029;
//            occupancy = 0.170293;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "CD1") {
            occupancy = 0.067;
//            occupancy = 0.258843;
        } else if (*neighboringAtom == "CD2") {
            occupancy = 0.076;
//            occupancy = 0.27568;
        } else if (*neighboringAtom == "CE1") {
            occupancy = 0.099;
//            occupancy = 0.31464;
        } else if (*neighboringAtom == "CE2") {
            occupancy = 0.048;
//            occupancy = 0.21909;
        } else if (*neighboringAtom == "CZ") {
            occupancy = 0.081;
//            occupancy = 0.284605;
        }
    } else if (*neighboringResi == "PRO") {
        if (*neighboringAtom == "N") {
            occupancy = 0.016;
//            occupancy = 0.12649;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.049;
//            occupancy = 0.221359;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.035;
//            occupancy = 0.18708;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.488;
//            occupancy = 0.698569;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.132;
//            occupancy = 0.363318;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.164;
//            occupancy = 0.404969;
        } else if (*neighboringAtom == "CD") {
            occupancy = 0.146;
//            occupancy = 0.382099;
        }
    } else if (*neighboringResi == "SER") {
        if (*neighboringAtom == "N") {
            occupancy = 0.245;
//            occupancy = 0.494975;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.022;
//            occupancy = 0.14832;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.072;
//            occupancy = 0.268328;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.150;
//            occupancy = 0.387298;
        } else if (*neighboringAtom == "OG") {
            occupancy = 0.568;
//            occupancy = 0.753657;
        }
    } else if (*neighboringResi == "THR") {
        if (*neighboringAtom == "N") {
            occupancy = 0.292;
//            occupancy = 0.54037;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.021;
//            occupancy = 0.1449137;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.601;
//            occupancy = 0.77524;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.054;
//            occupancy = 0.23238;
        } else if (*neighboringAtom == "OG1") {
            occupancy = 0.536;
//            occupancy = 0.73212;
        } else if (*neighboringAtom == "CG2") {
            occupancy = 0.193;
//            occupancy = 0.439317;
        }
    } else if (*neighboringResi == "TRP") {
        if (*neighboringAtom == "N") {
            occupancy = 0.429;
//            occupancy = 0.65498;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.046;
//            occupancy = 0.214476;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.031;
//            occupancy = 0.176068;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.569;
//            occupancy = 0.75432;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.049;
//            occupancy = 0.221359;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "CD1") {
            occupancy = 0.149;
//            occupancy = 0.386;
        } else if (*neighboringAtom == "CD2") {
            occupancy = 0.067;
//            occupancy = 0.25884;
        } else if (*neighboringAtom == "NE1") {
            occupancy = 0.349;
//            occupancy = 0.590762;
        } else if (*neighboringAtom == "CE2") {
            occupancy = 0.100;
//            occupancy = 0.31622;
        } else if (*neighboringAtom == "CE3") {
            occupancy = 0.027;
//            occupancy = 0.164316;
        } else if (*neighboringAtom == "CZ2") {
            occupancy = 0.072;
        } else if (*neighboringAtom == "CZ3") {
            occupancy = 0.071;
        } else if (*neighboringAtom == "CH2") {
            occupancy = 0.115;
        }
    } else if (*neighboringResi == "TYR") {
        if (*neighboringAtom == "N") {
            occupancy = 0.374;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.006;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.074;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.551;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.069;
        } else if (*neighboringAtom == "CG") {
            occupancy = 0.022;
        } else if (*neighboringAtom == "CD1") {
            occupancy = 0.083;
        } else if (*neighboringAtom == "CD2") {
            occupancy = 0.051;
        } else if (*neighboringAtom == "CE1") {
            occupancy = 0.068;
        } else if (*neighboringAtom == "CE2") {
            occupancy = 0.050;
        } else if (*neighboringAtom == "CZ") {
            occupancy = 0.036;
        } else if (*neighboringAtom == "OH") {
            occupancy = 0.635;
        }
    } else if (*neighboringResi == "VAL") {
        if (*neighboringAtom == "N") {
            occupancy = 0.482;
        } else if (*neighboringAtom == "CA") {
            occupancy = 0.013;
        } else if (*neighboringAtom == "C") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "O") {
            occupancy = 0.534;
        } else if (*neighboringAtom == "CB") {
            occupancy = 0.000;
        } else if (*neighboringAtom == "CG1") {
            occupancy = 0.085;
        } else if (*neighboringAtom == "CG2") {
            occupancy = 0.106;
        }
        // assign default occupancy for 1 for all other residues
    } else {
        occupancy = 1.0;
    };


    //cout << "Neighbor :" << *neighboringAtom << "|" << *neighboringResi << "|" << occupancy << endl;

    return occupancy;
}


unsigned int convertToAtomicNumber (const std::string atomType, const std::string resiName) {
    // Atomic Symbol is the first element of the name (right justified) " C  " vs "CA " vs " CA "
    //
    // Could have some problems with alpha-Carbon(C) and Calcium(Ca)
    //
    // Next element is the greek letter ( CB )
    // Next element is branch indicator ( CB1)
    // CB1 (carbon, beta, 1)
    // oxygen   atomType[0] =~ /\s/, atomType[1] == "O", atomType[2] =~ /[1-9PABGDEZH]/
    // nitrogen atomType[0] =~ /\s/, atomType[1] == "N", atomType[2] =~ /[1-9PABGDEZH]/
    // carbon   atomType[0] =~ /\s/, atomType[1] == "C", atomType[2] =~ /ABGDEZH/
    // calcium   atomType[0] == "C", atomType[1] == "A", atomType ==/^CA/, atomType[2] = /[1-9]/

    unsigned int atomicNumber = 6;

    string branch(1,atomType[3]); // get last element, cast const char* -> string

    if (isspace(atomType[0]) && (atomType[1] == 'C') ) {
        atomicNumber = 6;
    } else if (isspace(atomType[0]) && (atomType[1] == 'N')  ) {
        atomicNumber = 7;
    } else if (isspace(atomType[0]) && (atomType[1] == 'O') || ( resiName == "HOH" && atomType[1] == 'O')) {
        atomicNumber = 8;
    } else if (isspace(atomType[0]) && (atomType[1] == 'P') && boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 15;
    } else if (isspace(atomType[0]) && (atomType[1] == 'S') && boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 16;
    } else if (isspace(atomType[0]) && (atomType[1] == 'H') && boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 1;
    } else if (!isspace(atomType[0]) && (atomType[0] == 'H') && boost::regex_match(branch, columnNumBranch) ) {
        atomicNumber = 1;
    } else if (!isspace(atomType[0]) && (atomType[0] == 'H') && (boost::regex_match(string (1,atomType[1]), columnNumBranch) || (atomType[1] == 'O') )) {
        atomicNumber = 1;
    } else if (isspace(atomType[0]) && (atomType[1] == 'H') && (boost::regex_match(string (1,atomType[2]), columnNumBranch) || (atomType[2] == 'O') )) {
        atomicNumber = 1;
    } else if ((atomType[0] == 'H') && (atomType[1] == 'E') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 2;
    } else if ((atomType[0] == 'L') && (atomType[1] == 'I') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 3;
    } else if ((atomType[0] == 'B') && (atomType[1] == 'E') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 4;
    } else if (isspace(atomType[0]) && (atomType[1] == 'B') && boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 5;
    } else if (isspace(atomType[0]) && (atomType[1] == 'F') && boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 9;
    } else if ((atomType[0] == 'N') && (atomType[1] == 'E') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 10;
    } else if ((atomType[0] == 'N') && (atomType[1] == 'A') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 11;
    } else if ((atomType[0] == 'M') && (atomType[1] == 'G') &&  boost::regex_match(branch, columnAtomBranch) ) { //Mg
        atomicNumber = 12;
    } else if ((atomType[0] == 'A') && (atomType[1] == 'L') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 13;
    } else if ((atomType[0] == 'S') && (atomType[1] == 'I') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 14;
    } else if ((atomType[0] == 'C') && (atomType[1] == 'L') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 17;
    } else if ((atomType[0] == 'A') && (atomType[1] == 'R') &&  boost::regex_match(branch, columnAtomBranch) ) {
        atomicNumber = 18;
    } else if (isspace(atomType[0]) && (atomType[1] == 'K') && boost::regex_match(branch, columnAtomBranch) ) {	// K
        atomicNumber = 19;
    } else if ((atomType[0] == 'C') && (atomType[1] == 'A') &&  boost::regex_match(branch, columnAtomBranch) ) { // Ca
        atomicNumber = 20;
    } else if ((atomType[0] == 'M') && (atomType[1] == 'N') &&  boost::regex_match(branch, columnAtomBranch) ) { // Mn
        atomicNumber = 25;
    } else if ((atomType[0] == 'F') && (atomType[1] == 'E') &&  boost::regex_match(branch, columnAtomBranch) ) { // Fe
        atomicNumber = 26;
    } else if ((atomType[0] == 'C') && (atomType[1] == 'O') &&  boost::regex_match(branch, columnAtomBranch) ) { // Co
        atomicNumber = 27;
    } else if ((atomType[0] == 'N') && (atomType[1] == 'I') &&  boost::regex_match(branch, columnAtomBranch) ) { // Ni
        atomicNumber = 28;
    } else if ((atomType[0] == 'C') && (atomType[1] == 'U') &&  boost::regex_match(branch, columnAtomBranch) ) { // Cu
        atomicNumber = 29;
    } else if ((atomType[0] == 'Z') && (atomType[1] == 'N') &&  boost::regex_match(branch, columnAtomBranch) ) { // Zn
        atomicNumber = 30;
    } else if ((atomType[0] == 'S') && (atomType[1] == 'E')  ) { // SE
        atomicNumber = 34;
    } else if (isspace(atomType[0]) && (atomType[1] == 'O') && boost::regex_match(branch, columnAtomBranch) && (boost::regex_match(resiName, isWater)) ) { //HOH water
        atomicNumber = 99;
    } else {
        auto tempString = std::string(atomType);
        boost::algorithm::trim(tempString);
        if (tempString.length() == 1) {
            if (tempString[0] == 'N'){
                atomicNumber=7;
            } else if (tempString[0] == 'C'  ){
                atomicNumber=6;
            } else if ((tempString[0] == 'O')){
                atomicNumber=8;
            } else if (tempString[0] == 'P'){
                atomicNumber=15;
            } else {
                std::cout << "DEBUG: Couldn't find atomic number for atomType : " << atomType <<": " << atomType.length() << " | " << tempString << std::endl;
            }
        } else {
            std::string second(1,tempString[1]);
            boost::algorithm::trim(second);
            if (atomType[0] == 'N' && boost::regex_match(second, columnNumBranch) ){
                cout << "\n INCORRECT FORMAT => Guessing Atom => NITROGEN(7)   C1=> " << atomType[0] << " C2=> " << atomType[1] << " C3=> " << atomType[2] << " C4=> " << atomType[3] << endl;
                atomicNumber =7;
            } else if (atomType[0] == 'C' && boost::regex_match(second, columnNumBranch) ){
                cout << "\n INCORRECT FORMAT => Guessing Atom => CARBON(6)     C1=> " << atomType[0] << " C2=> " << atomType[1] << " C3=> " << atomType[2] << " C4=> " << atomType[3] << endl;
                atomicNumber =6;
            } else if ((atomType[0] == 'O' && boost::regex_match(second, columnNumBranch)) || (atomType[0] == 'O' && atomType[1] == 'P')){
                cout << "\n INCORRECT FORMAT => Guessing Atom => OXYGEN(8)     C1=> " << atomType[0] << " C2=> " << atomType[1] << " C3=> " << atomType[2] << " C4=> " << atomType[3] << endl;
                atomicNumber =8;
            } else if (atomType[0] == 'P'){
                cout << "\n INCORRECT FORMAT => Guessing Atom => PHOSPHATE(15) C1=> " << atomType[0] << " C2=> " << atomType[1] << " C3=> " << atomType[2] << " C4=> " << atomType[3] << endl;
                atomicNumber =15;
            } else {
                std::cout << "DEBUG: Couldn't find atomic number for atomType : " << atomType <<": " << atomType.length() << " | " << tempString << std::endl;
                if (tempString == "CA"){
                    atomicNumber = 6;
                    std::cout << "DEBUG : " << atomType << " => " << atomicNumber << std::endl;
                }
            }
        }
    }

    return atomicNumber;
}

float getAtomicMass(unsigned int atomicNumber){
    static std::map<unsigned int, float> masses = {
            {1, 1.0018}, {2, 4.0026}, {3, 6.9410}, {4, 9.0122}, {5, 10.811}, {6, 12.011},
            {7, 14.007}, {8, 15.999}, {9, 18.998}, {10, 20.180}, {11, 22.990}, {12, 24.305},
            {13, 26.982}, {14, 28.086}, {15, 30.974}, {16, 32.06}, {17, 35.45}, {18, 39.948},
            {19, 39.098}, {20, 40.078}, {21, 44.956}, {22, 47.867}, {23, 50.942}, {24, 51.996},
            {25, 54.938}, {26, 55.845}, {27, 58.933}, {28, 58.693}, {29, 63.546}, {30, 65.38},
            {31, 69.720}, {32, 72.590}, {33, 74.9216}, {34, 78.960}, {35, 79.904}, {36, 83.80},
            {37, 85.4678}, {38, 87.620}, {39, 88.9059}, {40, 91.22}, {41, 92.9064}, {42, 95.94},
            {43, 98}, {44, 101.070}, {45, 102.9055}, {46, 106.4}, {47, 107.868}, {48, 112.41},
            {49, 114.82}, {50, 118.69}, {51, 121.75}, {52, 127.60}, {53, 126.9045}, {54, 131.30},
            {55, 132.9054}, {56, 137.33}, {57, 138.9055}, {58, 140.120}, {59, 140.9077}, {60, 144.24},
            {61, 145}, {62, 150.4}, {63, 151.96}, {64, 157.25}, {65, 158.9254}, {66, 162.50},
            {67, 164.9304}, {68, 167.26}, {69, 168.9342}, {70, 173.04}, {71, 174.967}, {72, 178.49},
            {73, 180.9479}, {74, 183.850}, {75, 186.207}, {76, 190.20}, {77, 192.22}, {78, 195.090},
            {79, 196.9665}, {80, 200.590}, {81, 204.370}, {82, 207.20}, {83, 208.9804}, {84, 209},
            {85, 210}, {86, 222}, {87, 223}, {88, 226.0254}, {89, 227.0278}, {90, 232.0381},
            {91, 231.0359}, {92, 238.029}, {93, 237.0482}, {99, 18.01528}};

            auto it = masses.find(atomicNumber);

if (it == masses.end()){
    std::cout << " atomicNumber " << atomicNumber << std::endl;
}
            return (it != masses.end()) ? it->second : 0;
}

/**
	atomic scattering factor calculator using 5 gaussian appoximation. Acta Cryst (1995). A51, 416-431. Table 1.

	@param[in]	atomicNumber
	@param[in]	q
	@return     float atomic scattering form factor

 **/

float asf ( int atomicNumber, float q) {
    /*
	 * Taken from Waasmaier, D., and Kirfel, A. "New Analytical Scattering-Factor Functions for Free Atoms and Ions" Acta Cryst. (1995) A51, 416-431
	   The last element of the array contains the four gaussian terms for water taken from: Hajdu, F. Acta Cryst (1972). A28, 250. 	 */
    float f_q;
    static float coeffs[][11] =
            {   //a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5
                    {},
                    {0.413048,   0.294953,   0.187491,   0.080701,   0.023736,   0.000049,  15.569946,  32.398468,   5.711404,  61.889874,   1.334118}, //H
                    {0.732354,   0.753896,   0.283819,   0.190003,   0.039139,   0.000487,  11.553918,   4.595831,   1.546299,  26.463964,   0.377523}, //He
                    {0.974637,   0.158472,   0.811855,   0.262416,   0.790108,   0.002542,   4.334946,   0.342451,  97.102966, 201.363831,   1.409234}, //Li
                    {1.533712,   0.638283,   0.601052,   0.106139,   1.118414,   0.002511,  42.662079,   0.595420,  99.106499,   0.151340,   1.843093}, //Be
                    {2.085185,   1.064580,   1.062788,   0.140515,   0.641784,   0.003823,  23.494068,   1.137894,  61.238976,   0.114886,   0.399036}, //B
                    {2.657506,   1.078079,   1.490909,  -4.241070,   0.713791,   4.297983,  14.780758,   0.776775,  42.086842,  -0.000294,   0.239535}, //C
                    {11.893780,   3.277479,   1.858092,   0.858927,   0.912985,  -11.80490,   0.000158,  10.232723,  30.344690,   0.656065,   0.217287},//N
                    {2.960427,   2.508818,   0.637853,   0.722838,   1.142756,   0.027014,  14.182259,   5.936858,   0.112726,  34.958481,   0.390240}, //O
                    {3.511943,   2.772244,   0.678385,   0.915159,   1.089261,   0.032557,  10.687859,   4.380466,   0.093982,  27.255203,   0.313066}, //F
                    {4.183749,   2.905726,   0.520513,   1.135641,   1.228065,   0.025576,   8.175457,   3.252536,   0.063295,  21.813910,   0.224952}, //Ne
                    {4.910127,   3.081783,   1.262067,   1.098938,   0.560991,   0.079712,   3.281434,   9.119178,   0.102763, 132.013947,   0.405878}, //Na
                    {4.708971,   1.194814,   1.558157,   1.170413,   3.239403,   0.126842,   4.875207, 108.506081,   0.111516,  48.292408,   1.928171}, //Mg
                    {4.730796,   2.313951,   1.541980,   1.117564,   3.154754,   0.139509,   3.628931,  43.051167,   0.095960, 108.932388,   1.555918}, //Al
                    {5.275329,   3.191038,   1.511514,   1.356849,   2.519114,   0.145073,   2.631338,  33.730728,   0.081119,  86.288643,   1.170087}, //Si
                    {1.950541,   4.146930,   1.494560,   1.522042,   5.729711,   0.155233,   0.908139,  27.044952,   0.071280,  67.520187,   1.981173}, //P
                    {6.372157,   5.154568,   1.473732,   1.635073,   1.209372,   0.154722,   1.514347,  22.092527,   0.061373,  55.445175,   0.646925}, //S
                    {1.446071,   6.870609,   6.151801,   1.750347,   0.634168,   0.146773,   0.052357,   1.193165,  18.343416,  46.398396,   0.401005},
                    {7.188004,   6.638454,   0.454180,   1.929593,   1.523654,   0.265954,   0.956221,  15.339877,  15.339862,  39.043823,   0.062409},
                    {8.163991,   7.146945,   1.070140,   0.877316,   1.486434,   0.253614,  12.816323,   0.808945, 210.327011,  39.597652,   0.052821},
                    {8.593655,   1.477324,   1.436254,   1.182839,   7.113258,   0.196255,  10.460644,   0.041891,  81.390381, 169.847839,   0.688098},
                    {1.476566,   1.487278,   1.600187,   9.177463,   7.099750,   0.157765,  53.131023,   0.035325, 137.319489,   9.098031,   0.602102},
                    {9.818524,   1.522646,   1.703101,   1.768774,   7.082555,   0.102473,   8.001879,   0.029763,  39.885422, 120.157997,   0.532405},
                    {10.473575,   1.547881,   1.986381,   1.865616,   7.056250,   0.067744,   7.081940,   0.026040,  31.909672, 108.022842,   0.474882},
                    {11.007069,   1.555477,   2.985293,   1.347855,   7.034779,   0.065510,   6.366281,   0.023987,  23.244839, 105.774498,   0.429369},
                    {11.709542,   1.733414,   2.673141,   2.023368,   7.003180,  -0.147293,   5.597120,   0.017800,  21.788420,  89.517914,   0.383054},
                    {12.311098,   1.876623,   3.066177,   2.070451,   6.975185,  -0.304931,   5.009415,   0.014461,  18.743040,  82.767876,   0.346506},
                    {12.914510,   2.481908,   3.466894,   2.106351,   6.960892,  -0.936572,   4.507138,   0.009126,  16.438129,  76.987320,   0.314418},
                    {13.521865,   6.947285,   3.866028,   2.135900,   4.284731,  -2.762697,   4.077277,   0.286763,  14.622634,  71.966080,   0.004437},
                    {14.014192,   4.784577,   5.056806,   1.457971,   6.932996,  -3.254477,   3.738280,   0.003744,  13.034982,  72.554794,   0.265666},
                    {14.741002,   6.907748,   4.642337,   2.191766,  38.424042, -36.915829,   3.388232,   0.243315,  11.903689,  63.312130,   0.000397},
                    {15.758946,   6.841123,   4.121016,   2.714681,   2.395246,  -0.847395,   3.121754,   0.226057,  12.482196,  66.203621,   0.007238},
                    {16.540613,   1.567900,   3.727829,   3.345098,   6.785079,   0.018726,   2.866618,   0.012198,  13.432163,  58.866047,   0.210974},
                    {17.025642,   4.503441,   3.715904,   3.937200,   6.790175,  -2.984117,   2.597739,   0.003012,  14.272119,  50.437996,   0.193015},
                    {17.354071,   4.653248,   4.259489,   4.136455,   6.749163,  -3.160982,   2.349787,   0.002550,  15.579460,  45.181202,   0.177432},
                    {17.550570,   5.411882,   3.937180,   3.880645,   6.707793,  -2.492088,   2.119226,  16.557184,   0.002481,  42.164009,   0.162121},
                    {17.655279,   6.848105,   4.171004,   3.446760,   6.685200,  -2.810592,   1.908231,  16.606236,   0.001598,  39.917473,   0.146896},
                    {8.123134,   2.138042,   6.761702,   1.156051,  17.679546,   1.139548,  15.142385,  33.542667,   0.129372, 224.132507,   1.713368},
                    {17.730219,   9.795867,   6.099763,   2.620025,   0.600053,   1.140251,   1.563060,  14.310868,   0.120574, 135.771317,   0.120574},
                    {17.792040,  10.253252,   5.714949,   3.170516,   0.918251,   1.131787,   1.429691,  13.132816,   0.112173, 108.197029,   0.112173},
                    {17.859772,  10.911038,   5.821115,   3.512513,   0.746965,   1.124859,   1.310692,  12.319285,   0.104353,  91.777542,   0.104353},
                    {17.958399,  12.063054,   5.007015,   3.287667,   1.531019,   1.123452,   1.211590,  12.246687,   0.098615,  75.011948,   0.098615},
                    {6.236218,  17.987711,  12.973127,   3.451426,   0.210899,   1.108770,   0.090780,   1.108310,  11.468720,  66.684151,   0.090780},
                    {17.840963,   3.428236,   1.373012,  12.947364,   6.335469,   1.074784,   1.005729,  41.901382, 119.320541,   9.781542,   0.083391},
                    {6.271624,  17.906738,  14.123269,   3.746008,   0.908235,   1.043992,   0.077040,   0.928222,   9.555345,  35.860680, 123.552246},
                    {6.216648,  17.919739,   3.854252,   0.840326,  15.173498,   0.995452,   0.070789,   0.856121,  33.889484, 121.686691,   9.029517},
                    {6.121511,   4.784063,  16.631683,   4.318258,  13.246773,   0.883099,   0.062549,   0.784031,   8.751391,  34.489983,   0.784031},
                    {6.073874,  17.155437,   4.173344,   0.852238,  17.988686,   0.756603,   0.055333,   7.896512,  28.443739, 110.376106,   0.716809},
                    {6.080986,  18.019468,   4.018197,   1.303510,  17.974669,   0.603504,   0.048990,   7.273646,  29.119284,  95.831207,   0.661231},
                    {6.196477,  18.816183,   4.050479,   1.638929,  17.962912,   0.333097,   0.042072,   6.695665,  31.009790, 103.284348,   0.610714},
                    {19.325171,   6.281571,   4.498866,   1.856934,  17.917318,   0.119024,   6.118104,   0.036915,  32.529045,  95.037186,   0.565651},
                    {5.394956,   6.549570,  19.650681,   1.827820,  17.867832,  -0.290506,  33.326523,   0.030974,   5.564929,  87.130966,   0.523992},
                    {6.660302,   6.940756,  19.847015,   1.557175,  17.802427,  -0.806668,  33.031654,   0.025750,   5.065547,  84.101616,   0.487660},
                    {19.884502,   6.736593,   8.110516,   1.170953,  17.548716,  -0.448811,   4.628591,   0.027754,  31.849096,  84.406387,   0.463550},
                    {19.978920,  11.774945,   9.332182,   1.244749,  17.737501,  -6.065902,   4.143356,   0.010142,  28.796200,  75.280685,   0.413616},
                    {17.418674,   8.314444,  10.323193,   1.383834,  19.876251,  -2.322802,   0.399828,   0.016872,  25.605827, 233.339676,   3.826915},
                    {19.747343,  17.368477,  10.465718,   2.592602,  11.003653,  -5.183497,   3.481823,   0.371224,  21.226641, 173.834274,   0.010719},
                    {19.966019,  27.329655,  11.018425,   3.086696,  17.335455, -21.745489,   3.197408,   0.003446,  19.955492, 141.381973,   0.341817},
                    {17.355122,  43.988499,  20.546650,   3.130670,  11.353665, -38.386017,   0.328369,   0.002047,   3.088196, 134.907654,  18.832960},
                    {21.551311,  17.161730,  11.903859,   2.679103,   9.564197,  -3.871068,   2.995675,   0.312491,  17.716705, 152.192825,   0.010468},
                    {17.331244,  62.783924,  12.160097,   2.663483,  22.239950, -57.189842,   0.300269,   0.001320,  17.026001, 148.748993,   2.910268},
                    {17.286388,  51.560162,  12.478557,   2.675515,  22.960947, -45.973682,   0.286620,   0.001550,  16.223755, 143.984512,   2.796480},
                    {23.700363,  23.072214,  12.777782,   2.684217,  17.204367, -17.452166,   2.689539,   0.003491,  15.495437, 139.862473,   0.274536},
                    {17.186195,  37.156837,  13.103387,   2.707246,  24.419271, -31.586687,   0.261678,   0.001995,  14.787360, 134.816299,   2.581883},
                    {24.898117,  17.104952,  13.222581,   3.266152,  48.995213, -43.505684,   2.435028,   0.246961,  13.996325, 110.863091,   0.001383},
                    {25.910013,  32.344139,  13.765117,   2.751404,  17.064405, -26.851971,   2.373912,   0.002034,  13.481969, 125.836510,   0.236916},
                    {26.671785,  88.687576,  14.065445,   2.768497,  17.067781, -83.279831,   2.282593,   0.000665,  12.920230, 121.937187,   0.225531},
                    {27.150190,  16.999819,  14.059334,   3.386979,  46.546471, -41.165253,   2.169660,   0.215414,  12.213148, 100.506783,   0.001211},
                    {28.174887,  82.493271,  14.624002,   2.802756,  17.018515, -77.135223,   2.120995,   0.000640,  11.915256, 114.529938,   0.207519},
                    {28.925894,  76.173798,  14.904704,   2.814812,  16.998117, -70.839813,   2.046203,   0.000656,  11.465375, 111.411980,   0.199376},
                    {29.676760,  65.624069,  15.160854,   2.830288,  16.997850, -60.313812,   1.977630,   0.000720,  11.044622, 108.139153,   0.192110},
                    {30.122866,  15.099346,  56.314899,   3.540980,  16.943729, -51.049416,   1.883090,  10.342764,   0.000780,  89.559250,   0.183849},
                    {30.617033,  15.145351,  54.933548,   4.096253,  16.896156, -49.719837,   1.795613,   9.934469,   0.000739,  76.189705,   0.175914},
                    {31.066359,  15.341823,  49.278297,   4.577665,  16.828321, -44.119026,   1.708732,   9.618455,   0.000760,  66.346199,   0.168002},
                    {31.507900,  15.682498,  37.960129,   4.885509,  16.792112, -32.864574,   1.629485,   9.446448,   0.000898,  59.980675,   0.160798},
                    {31.888456,  16.117104,  42.390297,   5.211669,  16.767591, -37.412682,   1.549238,   9.233474,   0.000689,  54.516373,   0.152815},
                    {32.210297,  16.678440,  48.559906,   5.455839,  16.735533, -43.677956,   1.473531,   9.049695,   0.000519,  50.210201,   0.145771},
                    {32.004436,   1.975454,  17.070105,  15.939454,   5.990003,   4.018893,   1.353767,  81.014175,   0.128093,   7.661196,  26.659403},
                    {31.273891,  18.445440,  17.063745,   5.555933,   1.575270,   4.050394,   1.316992,   8.797154,   0.124741,  40.177994,   1.316997},
                    {16.777390,  19.317156,  32.979683,   5.595453,  10.576854,  -6.279078,   0.122737,   8.621570,   1.256902,  38.008820,   0.000601},
                    {16.839890,  20.023823,  28.428564,   5.881564,   4.714706,   4.076478,   0.115905,   8.256927,   1.195250,  39.247227,   1.195250},
                    {16.630795,  19.386616,  32.808571,   1.747191,   6.356862,   4.066939,   0.110704,   7.181401,   1.119730,  90.660263,  26.014978},
                    {16.419567,  32.738590,   6.530247,   2.342742,  19.916475,   4.049824,   0.105499,   1.055049,  25.025890,  80.906593,   6.664449},
                    {16.282274,  32.725136,   6.678302,   2.694750,  20.576559,   4.040914,   0.101180,   1.002287,  25.714146,  77.057549,   6.291882},
                    {16.289164,  32.807171,  21.095163,   2.505901,   7.254589,   4.046556,   0.098121,   0.966265,   6.046622,  76.598068,  28.096128},
                    {16.011461,  32.615547,   8.113899,   2.884082,  21.377867,   3.995684,   0.092639,   0.904416,  26.543257,  68.372963,   5.499512},
                    {16.070229,  32.641106,  21.489658,   2.299218,   9.480184,   4.020977,   0.090437,   0.876409,   5.239687,  69.188477,  27.632641},
                    {16.007385,  32.663830,  21.594351,   1.598497,  11.121192,   4.003472,   0.087031,   0.840187,   4.954467, 199.805801,  26.905106},
                    {32.563690,  21.396671,  11.298093,   2.834688,  15.914965,   3.981773,   0.801980,   4.590666,  22.758972, 160.404388,   0.083544},
                    {15.914053,  32.535042,  21.553976,  11.433394,   3.612409,   3.939212,   0.080511,   0.770669,   4.352206,  21.381622, 130.500748},
                    {15.784024,  32.454899,  21.849222,   4.239077,  11.736191,   3.922533,   0.077067,   0.735137,   4.097976, 109.464111,  20.512138},
                    {32.740208,  21.973675,  12.957398,   3.683832,  15.744058,   3.886066,   0.709545,   4.050881,  19.231543, 117.255005,   0.074040},
                    {15.679275,  32.824306,  13.660459,   3.687261,  22.279434,   3.854444,   0.071206,   0.681177,  18.236156, 112.500038,   3.930325},
                    {32.999901,  22.638077,  14.219973,   3.672950,  15.683245,   3.769391,   0.657086,   3.854918,  17.435474, 109.464485,   0.068033},
                    {33.281178,  23.148544,  15.153755,   3.031492,  15.704215,   3.664200,   0.634999,   3.856168,  16.849735, 121.292038,   0.064857},
                    {33.435162,  23.657259,  15.576339,   3.027023,  15.746100,   3.541160,   0.612785,   3.792942,  16.195778, 117.757004,   0.061755},
                    {15.804837,  33.480801,  24.150198,   3.655563,  15.499866,   3.390840,   0.058619,   0.590160,   3.674720, 100.736191,  15.408296},
                    {15.889072,  33.625286,  24.710381,   3.707139,  15.839268,   3.213169,   0.055503,   0.569571,   3.615472,  97.694786,  14.754303},
                    {33.794075,  25.467693,  16.048487,   3.657525,  16.008982,   3.005326,   0.550447,   3.581973,  14.357388,  96.064972,   0.052450}, //98, Einsteinium
                    {1.6607,      1.6277,      3.7734,    2.7903,     0,          0.1444,     0.3042,     5.1864,    12.7450,    30.7880,     0       }, //HOH, 99th element

            };

    q = q/(4.0*M_PI);  //need to convert scattering vector to 4PI*sin(theta)/lambda
    float * atom = coeffs[atomicNumber];
    static float q_squared = q*q;

    f_q =	  *(atom  )*expf(-*(atom+6 )*q_squared)
               + *(atom+1)*expf(-*(atom+7 )*q_squared)
               + *(atom+2)*expf(-*(atom+8 )*q_squared)
               + *(atom+3)*expf(-*(atom+9 )*q_squared)
               + *(atom+4)*expf(-*(atom+10)*q_squared)
               + *(atom+5);
    return f_q;
}



float kurtosis(float * residuals, int residualsSize){
//    float squared_x = 0.0f;
//    float fourth_x = 0.0f;
//    float third_x = 0.0f;
//    float x = 0.0f;
    float x_value;
    float divisor = 1.0f/(float)residualsSize;
    float k=0.0f;

    float meanit = 0.0f;

    for(int i=0; i<residualsSize; i++){
        x_value = residuals[i];
        meanit += x_value;
//        x += x_value;
//        x_squared = x_value*x_value;
//        squared_x += x_squared;
//        fourth_x += x_squared*x_squared;
//        third_x += x_squared*x_value;
    }


    meanit /= (float)residualsSize;
    float diff, diff2, top=0.0f, bottom=0.0f;
    for(int i=0; i<residualsSize; i++){
        x_value = residuals[i];
        diff = x_value - meanit;
        diff2 = diff*diff;
        top += diff2*diff2;
        bottom += diff2;
    }

    top /= (float)residualsSize;
    bottom /= (float)residualsSize;


//    float fourth_moment =0.0;
//    float divisor_2nd = divisor*divisor;
//    float divisor_4th = divisor_2nd*divisor_2nd;
//
//    fourth_moment = divisor*fourth_x + 4.0f*third_x*x*divisor_2nd - 4.0f*x*x*x*x*divisor_4th + 6.0f*squared_x*x*x*divisor*divisor_2nd + x*x*x*x*divisor_4th;
//
//    variance = divisor*squared_x - x*x*divisor_2nd;
//    k = fourth_moment/variance - 3.0f;
    k = top/(bottom*bottom) - 3;
    //cout << "kurtosis: " << k << " var " << variance << " fourth " << fourth_moment <<  endl;
    return k;
}


void calcFQ (const vector<float> &qvalues, set<int> atomList, float * f_q_array) {

    auto qvaluesSize = (int)qvalues.size();
    auto totalAtoms = (int)atomList.size();
    //
    // A linearized 2-d array => first atom in atomList
    for (int q = 0; q < qvaluesSize; q++) {
        // break out of iterator if atom already present
        int atomCount=0;
        for(auto & it : atomList){
            //std::cout << " Atom types " << it << std::endl;
            f_q_array[q*totalAtoms + atomCount] = asf(it, qvalues[q]);
            atomCount++;
        }
    }
}


float median(std::vector<float> * scores) {
    float median;
    size_t size = scores->size();

    sort(scores->begin(), scores->end());

    if (size  % 2 == 0) {
        median = ((*scores)[size / 2 - 1] + (*scores)[size / 2]) / 2;
    } else {
        median = (*scores)[size / 2];
    }

    return median;
}


void printInfo(std::string text) {
    std::cout << "  => " << text << std::endl;
}

void printError(std::string text) {
    std::cout << "*** ERROR => " << text << std::endl;
}


void logger(std::string description, std::string value) {
    unsigned int len = (description.size() > 40) ? 40 : description.size();
    printf("%40s :: %s\n", description.substr(0,len).c_str(), value.c_str());
}


bool validatePofRFile(std::string filename){

    bool hasBins = false;
    bool hasRealSpace = false;

    std::ifstream data (filename, std::ifstream::in);
    if (data.is_open()) {

        boost::regex format("REAL SPACE");
        boost::regex coefficient("BIN_[0-9]+", boost::regex::icase);
        boost::regex remarkFormat("REMARK");

        std::string line;

        while(!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if(boost::regex_search(line, format)){
                logger("READING PR VALUES", filename);
                hasRealSpace = true;
            }

            if(boost::regex_search(tempLine.at(0), remarkFormat) && boost::regex_search(line, coefficient) ){
                hasBins = true;
            }

            if (hasBins && hasRealSpace){
                break;
            }
        }
    }

    data.close();
    return (hasBins && hasRealSpace);
}

std::string formatNumber(unsigned int number) {
    char buffer [50];
    std::sprintf (buffer, "%d", number);
    return std::string(buffer);
}

std::string formatNumber(float number, int decimals) {
    char buffer [50];
    switch(decimals){
        case 1 :
            std::sprintf (buffer, "%.1f", number);
            break;
        case 2 :
            std::sprintf (buffer, "%.2f", number);
            break;
        case 3 :
            std::sprintf (buffer, "%.3f", number);
            break;
        default:
            std::sprintf (buffer, "%.4f", number);
    }

    return std::string(buffer);
}

void printAtomLine(FILE * pFile, unsigned int index, std::string chain, unsigned int residue_index, float x, float y, float z){
    std::string resid = std::to_string(residue_index);
    fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", index," CA ", "ALA", chain.c_str(), resid.c_str(), x, y, z );
}

DminType getDminValues(std::vector<vector3> & centered_coordinates){

    std::vector<float> dmin;
    const vector3 * pVec = centered_coordinates.data();
    unsigned int total_centered_coordinates = centered_coordinates.size();

    float dminSum = 0.0f;
    float dminSumSquare = 0.0f;
    for (unsigned int n=0; n < (total_centered_coordinates-1); n++) {
        const vector3 &pVec1 = *(pVec + n);
        float tempDmin = FLT_MAX;

        for (unsigned int m=0; m < n; m++) {
            const vector3 &pVec2 = *(pVec + m);
            float dis = (pVec1 - pVec2).length();
            if (dis < tempDmin) {
                tempDmin = dis;
            }
        }

        for (unsigned int m=n+1; m < total_centered_coordinates; m++) {
            const vector3 &pVec2 = *(pVec + m);
            float dis = (pVec1 - pVec2).length();
            if (dis < tempDmin) {
                tempDmin = dis;
            }
        }

        dmin.emplace_back(tempDmin);
        dminSum += tempDmin;
        dminSumSquare += tempDmin*tempDmin;
    }

    std::sort(dmin.begin(), dmin.end());

    DminType dmins;

    dmins.dmin_supremum = dmin[dmin.size()-1];
    dmins.dmin_infimum = dmin[0];
    dmins.stdev_dmin = std::sqrt(dminSumSquare/(float)dmin.size() - dminSum*dminSum/(float)(dmin.size()*dmin.size()));
    dmins.average_dmin = dminSum/(float)dmin.size();

    return dmins;
}