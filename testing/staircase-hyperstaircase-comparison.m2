-- Comparison of outputs between rPowerMinkowski, rPowerHR and rPowerStaircase
-- uses randomMonomialIdeal function from RandomIdeals package
-- V5 AUG 12, 2020 10:57
versionNumber := toString(5);

-- includes lego algorithm
-- used timing instead of elapsedTime

stairFilepath := "/Users/loser/Documents/Polymath/github/staircase/staircase-algo.m2";
hyperStairFilepath := "/Users/loser/Documents/Polymath/github/hyper-staircase/hyper-staircase.m2";

load stairFilepath;
load hyperStairFilepath;

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


printTestResults = (hyperStairTime, stairTime, outcome, i) -> (
  print "------------------------------------------------------------------------------";
  print concatenate("test ", toString(i), " of ", toString(testSize));
  print I;
  print concatenate("r = ", toString(r));
  print "\n";

  print "Hyper-staircase time:";
  print hyperStairTime;
  print "\n";

  print "Staircase time:";
  print stairTime;
  print "\n";
    
  print(outcome);
);



--------------------------------------------------------------------------------
print concatenate("\nRing vars: ", toString(flatten entries vars R), "\n");

for i from 1 to testSize do(
  -- step that randomizes r, if desired
  randomizeRationalPower();

  -- step that randomizes ideal, if desired. the function makeGeneratorDegreeList handles the exact randomness that user indicates
  if randomIdeal then (
    makeGeneratorDegreeList();
    I = randomMonomialIdeal(genDegList, R);
  );

  hyperStairTime := null;
  hyperStairAns := null;
  stairTime := null;
  stairAns := null;

  -- CALLING THE ALGORITHMS WE'RE TESTING, added randomness so that the order of calculation can be rotated, giving each algorithm a chacne to go first/second
  mySeed := random(1,30);

  if mySeed <= 15 then (
    hyperStairTime = timing hyperStairAns = rPowerHS(I,r);
    stairTime = timing stairAns = rPowerStaircase(I,r);
  ) else (
    stairTime = timing stairAns = rPowerStaircase(I,r);
    hyperStairTime = timing hyperStairAns = rPowerHS(I,r);
  );

  correctResult := hyperStairAns == stairAns;
  outcome := if correctResult then "Outcome: SAME RESULT\n" else "Outcome: DIFF RESULT, yikes\n";

  -- PRINTING THINGS ON TERMINAL
  printTestResults(hyperStairTime, stairTime, outcome, i);

  -- This would halt the program and warn user in the case of a failed test
  assert (correctResult);
);
