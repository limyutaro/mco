-- Comparison of outputs between rPowerMinkowski, rPowerHR and rPowerStaircase
-- uses randomMonomialIdeal function from RandomIdeals package
-- V5 AUG 12, 2020 10:57
versionNumber := toString(5);

-- includes lego algorithm
-- used timing instead of elapsedTime

stairFilepath := "/Users/loser/Documents/Polymath/m2-code/staircase-algo-nD-v8.m2";
hyperRectFilepath := "/Users/loser/Documents/Polymath/m2-code/hyperR-algo-v3.m2";
legoFilePath := "/Users/loser/Documents/Polymath/m2-code/lego-algo-v2.m2"
-- minkFilepath := "/Users/loser/Documents/Polymath/m2-code/minkowski/minkowskiAlg-w-locals.m2";

load stairFilepath;
load hyperRectFilepath;
load legoFilePath;
-- load minkFilepath;
loadPackage ("Polyhedra", Reload => true);
loadPackage ("RandomIdeals", Reload => true);
loadPackage ("SymbolicPowers", Reload => true);

-- specify number of tests to be conducted
testSize := 50;

-- specify ring R
R = QQ[w,x,y,z]; 

------------------ specify if necessary. If not, can ignore ------------------
-- if you want to fix the ideal, define it here:
I = ideal(x*y, y^2);

-- if you want to fix r, define it here:
r := 1;

-- if you want to fix the number of generators, define it here:
numOfGens := 8;

-- if you want to fix the degree of generators to be the same, define it here:
degOfAllGens := 8;

-- if you want to fix the degree and number of each generator, define it here:
-- (in which case, randumNumOfGens and randomDeg should both be false)
genDegList := {7, 5, 6, 7, 6};
--------------------------------------------------------------------------------
-- if you want the ideal tested to be in any way random, specify true below:
randomIdeal := true;

-- if you want random number of generators, be sure to specify the bounds:
randomNumOfGens := true;
maxNumOfGens := 8;
minNumOfGens := 1;

-- if you want random degrees on each generator, be sure to specify the bounds here:
randomDeg := true;
maxGenDegree := 10;
minGenDegree := 1;

-- if you want all degrees to be the same, put true below
fixDegsToBeSame := false;

-- specify limits for the "fineness" of the rational power, r.
--   randomR can be true of false, which would control whether you want r to be 
--     random or not. If true, testing will generate random rs
--   rLimits_0 is the max numerator size, rLimits_1 is the max denominator size
--   maxR is an upper bound for the r that the user is willing to have, for 
--     loop will help take care of the fact that random rational power 
--     generated is less than maxR
randomR := true;
rLimits := {100, 51};
maxR := 3;


-- file paths for test-cases, outputs, and data.

-- dimNumStr is just something I use to write them into files based on number of variables in a ring
-- be sure to indicate correctly in order to write test cases into the right file. Purpose is to segment test cases into smaller chunks
dimNumStr := toString(dim R);
testCaseFileNumber := toString();

testfilepath := concatenate("/Users/loser/Documents/Polymath/m2-code/testing/test-cases/", dimNumStr,"D/test-cases-", dimNumStr,"D-", testCaseFileNumber,".m2");
outputfilepath := concatenate("/Users/loser/Documents/Polymath/m2-code/testing/outputs/outputs-hr-stair-", dimNumStr, "D.txt");
datafilepath := concatenate("/Users/loser/Documents/Polymath/m2-code/testing/data-v", versionNumber, "/data-hr-stair-", dimNumStr, "D.txt");


--------------------------------------------------------------------------------


randomizeRationalPower = () -> (
  if randomR then (
    r = maxR + 1;
  
    -- ensures chosen r is does not exceed upper limit
    while r > maxR do (
      randomNume := random(1, rLimits_0);
      randomDenom := random(1, rLimits_1);
      r = randomNume / randomDenom;
    );
  );
);


--------------------------------------------------------------------------------


makeGeneratorDegreeList = () -> (
  -- randomizes the number of generators if needed
  if randomNumOfGens then (
    numOfGens = random(minNumOfGens, maxNumOfGens);
  );

  -- randomizes/fixes the degrees of each of the generators
  if randomDeg and not fixDegsToBeSame then (
    genDegList = {};

    for i from 1 to numOfGens do (
      genDegree := random(minGenDegree, maxGenDegree);
      genDegList = append(genDegList, genDegree);
    );
  ) else if fixDegsToBeSame then (
    -- if degree is random, then randomize a new degree
    if randomDeg then (
      degOfAllGens = random(minGenDegree, maxGenDegree);
    );

    -- initialize and empty list
    genDegList = {};

    -- make a list of the same number
    genDegList = toList(numOfGens:degOfAllGens);
  );
);


--------------------------------------------------------------------------------


printTestResults = (hyperRectTime, legoTime, stairTime, outcome, i) -> (
  print "------------------------------------------------------------------------------";
  print concatenate("test ", toString(i), " of ", toString(testSize));
  print I;
  print concatenate("r = ", toString(r));
  print "\n";

  print "Hyperrectangle time:";
  print hyperRectTime;
  print "\n";

  print "Lego time:";
  print legoTime;
  print "\n";

  print "Staircase time:";
  print stairTime;
  print "\n";

  -- print "Mink time:"
  -- print minkTime;
  -- print "\n";
    
  print(outcome);
);


--------------------------------------------------------------------------------


writeTestCases = (hyperRectAns) -> (
  testcasesfile := openOutAppend testfilepath;

  testcasesfile << concatenate("\nI = ", toString(I), ";") << endl;
  testcasesfile << concatenate("r = ", toString(r), ";") << endl;
  testcasesfile << concatenate("assert ( rPower(I, r) == ", toString(hyperRectAns), " );")<< endl;

  testcasesfile << close;
);


--------------------------------------------------------------------------------


writeOutputs = (hyperRectTime, legoTime, stairTime, outcome) -> (
  outputsfile := openOutAppend outputfilepath;

  outputsfile << "------------------------------------------------------------------------------\n" << endl;

  outputsfile << I << endl;
  outputsfile << concatenate("r = ", toString(r), "\n") << endl;

  outputsfile << "Hyperrectangle:" << endl;
  outputsfile << hyperRectTime << endl;
  outputsfile << endl;
  
  outputsfile << "Lego:" << endl;
  outputsfile << legoTime << endl;
  outputsfile << endl;

  outputsfile << "Staircase:" << endl;
  outputsfile << stairTime << endl;
  outputsfile << endl;

  outputsfile << outcome << endl;

  outputsfile << close;
);


--------------------------------------------------------------------------------


writeData = (hyperRectTime, legoTime, stairTime, hyperRectAns) -> (
  datafile := openOutAppend datafilepath;

  -- size of ring, r, ideal
  numberOfVariablesInRing := concatenate(dimNumStr, ",");
  rationalPower := concatenate(toString(r), ",");
  stringOfIdeal := concatenate(toString(I), ",");
    
  -- height and big height of ideal, not adjusted by rational power r
  idealHeight := concatenate(toString(codim I), ",");
  idealBigHeight := concatenate(toString(bigHeight I), ",");

  -- timing for hyperrectangle algorithm
  runtimeForHyperRect := toString(take(hyperRectTime, 1));
  runtimeForHyperRect = substring(runtimeForHyperRect, 5);
  runtimeForHyperRect = replace("}", "", runtimeForHyperRect);
  runtimeForHyperRect = concatenate(runtimeForHyperRect, ",");

  -- timing for lego algorithm
  runtimeForLegoRect := toString(take(legoTime, 1));
  runtimeForLegoRect = substring(runtimeForLegoRect, 5);
  runtimeForLegoRect = replace("}", "", runtimeForLegoRect);
  runtimeForLegoRect = concatenate(runtimeForLegoRect, ",");

  -- timing for staircase algorithm
  runtimeForStaircase := toString(take(stairTime, 1));
  runtimeForStaircase = substring(runtimeForStaircase, 5);
  runtimeForStaircase = replace("}", "", runtimeForStaircase);
  runtimeForStaircase = concatenate(runtimeForStaircase, ",");

  -------------- INTERMEDIATE DATA, NOT BEING WRITTEN ON DATA FILE -------------
  -- used a helper to get scaled exponents and min/max
  infoAboutGenerators := getExtremeValuesTest(I,r);
  
  -- reformat scaled exponents
  generatorsList := infoAboutGenerators_0;
  scaledGenList := apply(r * generatorsList, i -> apply(i, ceiling));
  
  -- get degree of each generator
  degreeOfEachGenerator = apply(scaledGenList, sum);
  
  -- extreme values and their differences
  minValues := infoAboutGenerators_1;
  maxValues := infoAboutGenerators_2;
  sideLengthsOfHyperRect := apply(maxValues, minValues, (i,j) -> i - j);
  ------------------------------------------------------------------------------

  -- describing the generators themselves. Thus, it is information regarding the vertices, not so much the hyperrectangle region of our search
  numberOfGenerators := concatenate(toString(numgens I), ",");
  
  largestDegree := concatenate(toString(max degreeOfEachGenerator), ",");
  smallestDegree := concatenate(toString(min degreeOfEachGenerator), ",");
  totalDegree := sum degreeOfEachGenerator;
  averageDegree := concatenate(toString(totalDegree / length degreeOfEachGenerator), ",");
  totalDegree = concatenate(toString(totalDegree), ",");

  numberOfVariablesInEachGenerator := apply(scaledGenList, i -> number(i, j -> j != 0));
  largestNumberOfVariablesInAGenerator := concatenate(toString(max numberOfVariablesInEachGenerator), ",");
  smallestNumberOfVariablesInAGenerator := concatenate(toString(min numberOfVariablesInEachGenerator), ",");
  averageNumberOfVariablesInEachGenerator := concatenate(toString(sum numberOfVariablesInEachGenerator / length numberOfVariablesInEachGenerator), ",");

  volumeOfNewtonPolytopeOfRPower := concatenate(toString(volume newtonPolytope sum flatten entries gens hyperRectAns), ",");

  -- describing the hyperrectangle
  volumeOfHyperRect := concatenate(toString(fold(sideLengthsOfHyperRect, (i,j) -> i * j)), ",");
  longestSideLength := concatenate(toString(max sideLengthsOfHyperRect), ",");
  shortestSideLength := concatenate(toString(min sideLengthsOfHyperRect), ",");
  totalSideLength := sum sideLengthsOfHyperRect;
  averageSideLength := concatenate(toString(totalSideLength / length sideLengthsOfHyperRect), ",");
  totalSideLength = concatenate(toString(totalSideLength), ",");
  
  -- concatenation of all the datum found earlier into comma seperated string
  datum := concatenate(numberOfVariablesInRing, runtimeForHyperRect, runtimeForLegoRect, runtimeForStaircase, rationalPower, numberOfGenerators, largestDegree, smallestDegree, totalDegree, averageDegree, largestNumberOfVariablesInAGenerator, smallestNumberOfVariablesInAGenerator, averageNumberOfVariablesInEachGenerator, volumeOfHyperRect, longestSideLength, shortestSideLength, totalSideLength, averageSideLength, idealHeight, idealBigHeight, stringOfIdeal, volumeOfNewtonPolytopeOfRPower);

  -- write data into file
  datafile << datum << endl << close;
);


--------------------------------------------------------------------------------


getExtremeValuesTest = (I, r) ->
(
  generatorArray := apply(apply(flatten entries gens I, exponents), flatten);

  -- INITIALIZE LIST OF START AND END POINTS TO A LIST OF SIZE DIMENSIONS, EACH ELEMENT BEING THE FIRST GENERATOR IN generatorArray
  minValues := r * generatorArray_0;
  maxValues := minValues;

  restOfGenerators := drop(generatorArray, 1);

  -- FOR EACH x_i, FIND THE POINTS WITH SMALLEST AND LARGEST x_i EXPONENTS
  scan(restOfGenerators, point -> 
      (
        newPoint := r * point;

        -- FOR EACH VARIABLE OF newPoint, CHECK IF IT IS A MAX OR MIN
        for i from 0 to (numgens(R) - 1) do(
          candidateValue := newPoint_i;

          if (candidateValue < minValues_i) then (
            minValues = replace(i, candidateValue, minValues);
          ) else if (candidateValue > maxValues_i) then (
            maxValues = replace(i, candidateValue, maxValues)
          ); 
        );
      )
    );

  -- Take ceiling for all values, at every coordinate
  minValues = apply(minValues, value -> ceiling(value));
  maxValues = apply(maxValues, value -> ceiling(value));
  
  return {generatorArray, minValues, maxValues}
);


--------------------------------------------------------------------------------
print concatenate("\nRing vars: ", toString(flatten entries vars R));

print "\n";

for i from 1 to testSize do(
  -- step that randomizes r, if desired
  randomizeRationalPower();

  -- step that randomizes ideal, if desired. the function makeGeneratorDegreeList handles the exact randomness that user indicates
  if randomIdeal then (
    makeGeneratorDegreeList();
    I = randomMonomialIdeal(genDegList, R);
  );

  hyperRectTime := null;
  hyperRectAns := null;
  legoTime := null;
  legoAns := null;
  stairTime := null;
  stairAns := null;
  -- minkTime := null;
  -- minkAns := null;

  -- CALLING THE ALGORITHMS WE'RE TESTING, added randomness so that the order of calculation can be rotated, giving each algorithm a chacne to go first/second
  mySeed := random(1,30);

  if mySeed <= 10 then (
    hyperRectTime = timing hyperRectAns = rPowerHR(I,r);
    stairTime = timing stairAns = rPowerStaircase(I,r);
    legoTime = timing legoAns = rPowerLego(I,r);
    -- minkTime = elapsedTiming minkAns = rPowerMinkowski(I,r);  
  ) else if mySeed <= 20 then (
    stairTime = timing stairAns = rPowerStaircase(I,r);
    legoTime = timing legoAns = rPowerLego(I,r);
    hyperRectTime = timing hyperRectAns = rPowerHR(I,r);
  ) else (
    legoTime = timing legoAns = rPowerLego(I,r);
    hyperRectTime = timing hyperRectAns = rPowerHR(I,r);
    stairTime = timing stairAns = rPowerStaircase(I,r);
  );

  correctResult := hyperRectAns == stairAns;
  correctResult = correctResult and (hyperRectAns == legoAns);
  outcome := if correctResult then "Outcome: SAME RESULT\n" else "Outcome: DIFF RESULT, yikes\n";

  -- PRINTING THINGS ON TERMINAL
  printTestResults(hyperRectTime, legoTime, stairTime, outcome, i);

  -- WRITING THINGS IN TO TEST CASES FILE
  writeTestCases(hyperRectAns);

  -- WRITING THINGS IN TO OUTPUTS FILE
  writeOutputs(hyperRectTime, legoTime, stairTime, outcome);

  -- WRITING THINGS IN TO DATA FILE
  writeData(hyperRectTime, legoTime, stairTime, hyperRectAns);

  -- This would halt the program and warn user in the case of a failed test
  assert (correctResult);
);
