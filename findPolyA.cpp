#include <Rcpp.h>
#include <fstream>
#include <sstream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

//////////////////////////////////////////////////////////////////////////////////////////
//  Help functions
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * This function calculates the sum of entries in an array between start and end idx
 * 
 * @param vec The vector to calculate the sum from
 * @param start The idx from where in vec to start
 * @param end The idx where to end in vec, should not be outOfBones
 * @return int The computed sum
 */
int sum(vector<int> vec, int start, int end){
  int sum = 0;
  for(int i = start; i <= end; i++) sum += vec[i];
  return sum;
}

/**
 * This function is a simple implementation of "Run Length Encoding", also known as "rle".
 * The function simply computes the lengths and values of runs of equal values in the
 * given input vector. 
 * 
 * The region of interest is defined by lowerBoundary and upperBoundary, in this case
 * the region is computed to be 1 otherwise 0. Regions cumputed to 0 can be found if
 * idx % 2 = 0 and regions to 1 can be found if idx % 2 = 1 of the returned vector.
 * 
 * @param eventMean The event mean vector to compute
 * @param lowerBoundary The minimum height of positive region
 * @param upperBoundary The maximum height of positive region
 * @return The computed rle vector containing the sum of idx with same value
 */
vector<int> rle(vector<float> eventMean, float lowerBoundary, float upperBoundary){
  vector<int> result; int sum = 1; bool state = false;
  
  for(int i = 0; i < eventMean.size(); i++){
    float value = eventMean[i];
    bool status = (value > lowerBoundary && value < upperBoundary);
    
    // Add dummy if first state is positive, such that first state always stays negative
    if(i == 0 && status){
      result.push_back(0);
      state = true;
    }
    
    if(status == state)
      sum++;
    else{
      state = !state;
      result.push_back(sum);
      sum = 1;
    }
  }
  
  return result;
}

/**
 * This function finds the beginning of a possible polyA region with high porbability on
 * most samples. The function is analyzing the mean using rle on event data.
 * 
 * @param eventMean The vector from event data containing all the mean values
 * @param lowerBoundary Used for rle, setting the lower boundary for mean in polyA region
 * @param upperBoundary Used for rle, setting the upper boundary for mean in polyA region
 * @param lowerRegionTimeThreshold Value for length of lower region before polyA
 * @param noisyBeginningThreshold Value for length of noisy region nearby polyA region
 * 
 * @return Index (based on eventMean) of the beginnng of a possible polyA region 
 */
int findBeginningOfPolyARegion(vector<float> eventMean, float lowerBoundary, 
          float upperBoundary, int lowerRegionTimeThreshold, int noisyBeginningThreshold){
  
  // PolyA region is assumed to be between lowerBoundary and upperBoundary
  vector<int> htable = rle(eventMean, lowerBoundary, upperBoundary);
  if(htable.empty()) return -1;
  
  // -- Print for DEBUG -- //
  // cout << "hTable:" << endl;
  // for(int i = 0; i < htable.size(); i++){
  //   cout << htable[i] << ", ";
  // }
  // cout << endl << endl;
  
  int possibleStartIdx = 0;
  for(int i = 2; i < htable.size()-5; i++){
    
    // Before polyA region, I assume that the event mean stays lower than lowerBoundary 
    // over a long time threshold
    if(i % 2 == 0 && htable[i] > lowerRegionTimeThreshold){

      // Check for noisy beginning (2 steps) with a threshold value.
      // This also assumes the minimum length of a part of polyA region.
      int n = noisyBeginningThreshold;
      if(htable[i+1] < n){
        if(htable[i+2] < n && (htable[i+3] > n || (htable[i+4] < n && htable[i+5] > n ))){
          possibleStartIdx = sum(htable, 0, i);
          break;
        }
        continue;
      }

      possibleStartIdx = sum(htable, 0, i);
      break;
    }
  }
  return possibleStartIdx;
}

/** 
 * This functions finds gaps in eventData with no moves.
 * 
 * These gaps are based on a threshold value. In the eventMove the function expects to 
 * find big parts with only zeros. If such a part is bigger than the threshold value, then 
 * this part becomes a gap.
 * 
 * Use start value from where to start, useful to avoid not necessary gap finding.
 * 
 * The function returns the idx's in event data where a gap begins and ends. 
 * 
 * Smaller threshold will result in more founded gaps.
 * 
 * @param eventMove The eventMove array from event data table.
 * @param start The idx from where in eventMove you want starting gap finding.
 * @param gapTreheshold The threshold value (minimum length of a zero part, to get a gap)
 * @return Vector containing all cut idx's, idx in event data
 */
vector< vector<int> > findGapsInEventMove(vector<int> eventMove, int start, int gapThreshold){
  vector< vector<int> > result; int cutter = 0; 
  for(int i = start; i < eventMove.size(); i++){
    // Beginning or middle of zero part
    if(eventMove[i] == 0){
      cutter++;
    }
    // End of zero part or something other
    else{
      if(cutter > 0 && cutter > gapThreshold){
        vector<int> tmp;
        tmp.push_back(i-cutter); tmp.push_back(i);
        result.push_back(tmp);
      }
      cutter = 0;
    }
  }
  return result;
}

/**
 * This function is filtering the gaps obtained by findGapsInEventMove() function. 
 * Using a moveThreshold the function is removing gaps with lower moves than the 
 * threshold value to the next gap. The function begins on lower regions and is stopping
 * at the gap containing the possible end of the polyA region, where the moves to next
 * gap becomes much bigger than the threshold value.
 * 
 * @param gaps Vector containing all the gaps obtained by findGapsInEventMove() function
 * @param eventMove The eventMove array from event data table
 * @param moveThreshold The minimum amount of moves between two gaps.
 */
int filterLowerRegions(vector< vector<int> > gaps, vector<int> eventMove, int moveThreshold){
  
  // Iterating throw gaps and checking if number of moves between gaps becomes higher
  // than some threshold
  int cutIdx = 0;
  for(int i = 0; i < gaps.size()-1; i++){

    int tmpE = gaps[i][1]; int tmpS = gaps[i+1][0]; int sum = 0;
  
    for(int j = tmpE; j < tmpS; j++) 
      if(eventMove[j] > 0) 
        sum++;
      
    if(sum < moveThreshold){
      // -- debug -- //
      //cout << sum << ",, "; 
      
      cutIdx = i+1;
      continue;
    }
    else {
      // -- debug -- //
      //cout << "M: " << sum << " From: " << tmpE*15 << " To: " << tmpS*15 << endl;
      
      break;
    }
  }
  return cutIdx;
}

/**
 * This function puts everything together and evaluating the result. The function returns
 * the beginning of the polyA region, the end and the length.
 * 
 * @param gaps Vector containing all the gaps obtained by findGapsInEventMove() function
 * @param eventMean The eventMean array from event data table
 * @param pStart Possible start of polyA region, from findBeginningOfPolyARegion() function
 * @param cutIdx Possible end region of polyA region, from filterLowerRegions() function
 * @param lowerBoundary The lower boundary for mean in polyA region
 * @param upperBoundary The upper boundary for mean in polyA region
 * 
 * @return Vector of length three, containing start, end and length of polyA as event Idx
 */
vector<int> evaluateGapsForPolyA(vector< vector<int> > gaps, vector<float> eventMean, 
                        int pStart, int cutIdx, float lowerBoundary, float upperBoundary){
  vector<int> result;
  
  // Caluclate means of mean for the gap
  vector<float> gapsMean;
  for(vector<int> vec : gaps){
    float sum = 0;
    for(int i = vec[0]; i < vec[1]; i++)
      sum += eventMean[i];
    gapsMean.push_back(sum/float(vec[1] - vec[0]));
  }
  
  // Check if actual gap at cutIdx is part of PolyA
  int upperIdx = -1, lowerIdx = -1;
  if(gapsMean[cutIdx] < lowerBoundary || gapsMean[cutIdx] > upperBoundary){
    
    // Check of cutIdx is too high and check if polyA region may be lower
    if(cutIdx > 0){
      for(int i = cutIdx - 1; i >= 0; i--){
        if(gapsMean[i] > lowerBoundary && gapsMean[i] < upperBoundary){
          lowerIdx = gaps[i][0];
          upperIdx = gaps[i][1];
          cutIdx = i;
          break;
        }
      }
    }
    
    // -- TODO: If not part of polyA, check alternative approach? -- //
    // --               (reject sample by now)                    -- //
    else{
      result.push_back(-1);
      result.push_back(-1); 
      result.push_back(-1);
      return result;
    }
    
  }
  else{
    lowerIdx = gaps[cutIdx][0];
    upperIdx = gaps[cutIdx][1];
  }
  
  // Evaluate drops in PolyA by checking upper gaps in neightborhood
  for(int i = cutIdx + 1; i < gaps.size(); i++){
    if(gapsMean[i] > lowerBoundary && gapsMean[i] < upperBoundary && (gaps[i][0] - upperIdx) < 100)
      upperIdx = gaps[i][1];
    else
      break;
  }
  
  // Evaluate drops in PolyA by checking lower gaps in neightborhood
  for(int i = cutIdx - 1; i >= 0; i--){
    if(gapsMean[i] > lowerBoundary && gapsMean[i] < upperBoundary && (lowerIdx - gaps[i][1]) < 100)
      lowerIdx = gaps[i][0];
    else
      break;
  }
  
  // Evaluate pStart vs lowerIdx, decide the true start of PolyA region
  if(pStart > 0 && pStart < lowerIdx && lowerIdx-pStart < 100){
    lowerIdx = pStart;
  }
  
  result.push_back(lowerIdx);
  result.push_back(upperIdx); 
  result.push_back(upperIdx-lowerIdx);
  
  // -- debug -- //
  //cout << "low: " << lowerIdx*15 << " high: " << upperIdx*15 << endl;
  
  return result;
}

/**
 * This function calculates the dwelltime for one basepair in the sample
 * 
 * @param eventMove The eventMove array from event data table
 * @param endPolyARegion The event idx where the polyA regions ends
 * @param tailEventStart The last value in eventStart + eventScope
 * @return float that is representing the dwelltime for one basepair in the sample
 */
float computeDwelltimePerBasepair(vector<int> eventMove, int endPolyARegion, int tailEventStart, int eventLength){
  float moves = sum(eventMove, endPolyARegion, eventMove.size()-1);
  return (float(tailEventStart - endPolyARegion*eventLength))/moves;
}

//////////////////////////////////////////////////////////////////////////////////////////
//  Main/Exported functions
//////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List findPolyA(List data) {
  
  //  ** All the treshold values ** //
  
  // I expect polyA region between 0.3 and 1.5 in event mean
  float lowerBoundary = 0.3;
  float upperBoundary = 1.5;
  
  // Before polyA, i expect a big lower region (less than 0.3) with length about 50 and
  // a noisy beginning up to 10 in length
  int lowerRegionTimeThreshold = 50;
  int noisyBeginningThreshold = 10;
  
  // The minimum length/size of a gap is assumed to be 25, also the moves between gaps
  // should be more than 25 moves to be part of tx (after polyA region)
  int gapThreshold = 25;
  int moveThreshold = 25; 
  
  
  // ** The main findPolyA algorithm ** //
  
  vector<int> sPolyAdw, ePolyAdw, lPolyAdw;
  vector<float> lPolyAbp, dwPerBasePair;
  
  for(List s : data){
    vector<int> eventStart = as< vector<int> >(s["eventStart"]);
    vector<int> eventMove = as< vector<int> >(s["eventMove"]);
    vector<float> eventMean = as< vector<float> >(s["eventMean"]);
    int eventLength = as <int>(s["eventLength"]);
    
    // Check if event data of sample begins with 0 dwelltime, reject sample otherwise
    if(eventStart[0] > 0){
      sPolyAdw.push_back(-2); ePolyAdw.push_back(-2); lPolyAdw.push_back(-2); 
      lPolyAbp.push_back(-2); dwPerBasePair.push_back(-2);
      continue; 
    }
    
    // Find idx in event data where the polyA region may begin, no beginning rejects sample
    int pStart = findBeginningOfPolyARegion(eventMean, lowerBoundary, upperBoundary, 
                                         lowerRegionTimeThreshold, noisyBeginningThreshold);
    if(pStart < 0){
      sPolyAdw.push_back(-3); ePolyAdw.push_back(-3); lPolyAdw.push_back(-3); 
      lPolyAbp.push_back(-3); dwPerBasePair.push_back(-3);
      continue; 
    }
    
    // -- debug -- //
    //cout << "pStart: " << pStart << endl;
    
    // Find gaps of move, no gaps results in rejecting the sample
    vector< vector<int> > gaps = findGapsInEventMove(eventMove, pStart, gapThreshold);
    if(gaps.empty()){
      sPolyAdw.push_back(-4); ePolyAdw.push_back(-4); lPolyAdw.push_back(-4); 
      lPolyAbp.push_back(-4); dwPerBasePair.push_back(-4);
      continue; 
    }
    
    // -- debug -- //
    // for(vector<int> i : gaps){
    //   for(int j : i){
    //     cout << j*15 << ", ";
    //   }
    //   cout << endl;
    // }
    
    // Filter moves and find a possible end idx
    int cutIdx = filterLowerRegions(gaps, eventMove, moveThreshold);
    
    // -- debug -- //
    //cout << "cutIdx: " << cutIdx << endl;
    
    // Evaluate everything, finding start, end and length in dwelltime
    vector<int> polyA = evaluateGapsForPolyA(gaps, eventMean, pStart, cutIdx, lowerBoundary, upperBoundary);
    if(polyA[0] < 0){
      sPolyAdw.push_back(-1); ePolyAdw.push_back(-1); lPolyAdw.push_back(-1); 
      lPolyAbp.push_back(-1); dwPerBasePair.push_back(-1);
      continue; 
    }
    
    // -- debug -- //
    // for(int i : result){
    //   cout << i << ", ";
    // }
    // cout << endl;
  
    // Calculate DwellTime per one basepair for sample
    int tailStart = eventStart[eventStart.size()-1] + eventLength;
    float dwPerBp = computeDwelltimePerBasepair(eventMove, polyA[1], tailStart, eventLength);
    
    // -- debug -- //
    //cout << "dwPerBp: " << dwPerBp << endl;
  
    // Find length of polyA in basepair for sample and prepare result
    sPolyAdw.push_back(polyA[0] * eventLength); 
    ePolyAdw.push_back(polyA[1] * eventLength); 
    lPolyAdw.push_back(polyA[2] * eventLength); 
    lPolyAbp.push_back(float(polyA[2]) * float(eventLength) / dwPerBp); 
    dwPerBasePair.push_back(dwPerBp);
    
  }

  return List::create(
    Named("startPolyA_dwelltime") = sPolyAdw,
    Named("endPolyA_dwelltime") = ePolyAdw,
    Named("lengthPolyA_dwelltime") = lPolyAdw,
    Named("lengthPolyA_BasePair") = lPolyAbp,
    Named("DwelltimePerBasepair") = dwPerBasePair
  );
}

// [[Rcpp::export]]
List analyzeFastQ(string path){
  
  vector<string> traceNames;
  vector<string> path2fast5;
  
  string line;
  ifstream fastq (path);
  if(fastq.is_open()){
    while(getline (fastq, line)){
      
      if(line[0] != '@') continue;
      
      line.erase(0, 1);
      istringstream lineStream(line);
      vector<string> split; string item;
      while(getline(lineStream, item, ' ') ){
        split.push_back(item);
      }
      
      traceNames.push_back(split[0]);
      path2fast5.push_back(split[2]);
      
    }
    fastq.close();
  }
  else{
    cout << "[ERROR]  Could not find fastq file!" << endl;
  }
  
  return List::create(
    Named("traceName") = traceNames,
    Named("fast5_path") = path2fast5
  );
}