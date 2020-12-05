-- n-D STAIRCASE ALGORITHM V7 w path prints AUG 5 14:46

-- failed at some cases on 4 and 5 dimensions
-- Change to account for when waterfall begins at a point far too much in interiror of NP, which leasd to certain generators not caught
-- added function moveStartPointToBoundary to make the starting point run to edge before performing algorithms
-- changed comments at waterfall move up/down explanation
-- added better prints for whether the extreme is min or max
-- added extra prints at the end to print the ideals out in both forms.


-- In the description of the functions below, we will be loose by interchanging 
--   coordinate with variable and value with exponent, because they are 
--   conceptually equivalent.


--------------------------------------------------------------------------------
-- Place to define the ring(s), ideal(s) and rational power(s)

-- R = QQ[w,x,y,z];
-- I = ideal(x * y^2, y * z, x^3 * z^3);
-- I = ideal(x * y, y * z, x * z);
-- I = ideal(x^2*y^2, x^3*z^3, y^3, y^2*z^5);
-- I = ideal(x^2*z, z^2*x, y^2*z^2);
-- I = ideal(x^3*y^2, x*z^2, x^4*z, y^2*z^3);
-- I = ideal(x^3*y^2, x*z^2, x^4*z, y*z^3, y^2*z);
-- I = ideal (x^7*y^3, x^5*z, y^5);
-- I = ideal(w^2, w*z, x^2*y^3);
-- I = ideal(w*y, w*x, x^4*z);
-- I = ideal(w*z, x*y, x^3*z^2 , x^2 *z^3 );
-- r = 2;


-- rPowerStaircase(I,r)

--------------------------------------------------------------------------------
-- Loading relevant package needed

loadPackage ("Polyhedra", Reload => true);

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function given by M2 intro lecture

-- Appears in getHalfSpaces function

-- input: an ideal I

-- output: the Newton Polyhedron of the ideal
--------------------------------------------------------------------------------

Npolyhedron = I ->
(
  C = posHull id_(ZZ^(dim ring I));

  NP = C + newtonPolytope sum flatten entries gens I;

  return NP
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that finds the extreme values of the staircase algorithm, giving us
--   boundary conditions needed for the algorithm to terminate.

-- Appears in rPowerStaircase function

-- input: an ideal I and a rational number r

-- output: a list containing three lists. 
--         The first list represents the start point of our algorithm, which --           can be any point with the lowest value at the first coordinate 
--           (smallest exponent at the x_1 variable, among all minimal 
--           generators of I).
--         The second is the list of min values, the ith element of the list is 
--           the min value at the ith coordinate, among all minimal generators --           of I. 
--         The third is the list of max values, the ith element of the list is 
--           the max value at the ith coordinate, among all minimal generators --           of I. 
--------------------------------------------------------------------------------

getExtremeValues = (I, r) ->
(
  generatorArray := apply(apply(flatten entries gens I, exponents), flatten);

  -- INITIALIZE LIST OF START AND END POINTS TO A LIST OF SIZE DIMENSIONS, EACH ELEMENT BEING THE FIRST GENERATOR IN generatorArray
  dummyGenerator := r * generatorArray_0;
  startPoint := dummyGenerator;
  minValues := {};

  for i from 0 to numgens(R) - 1 do(
    minValues = append(minValues, dummyGenerator_i)
  );

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
            if i == 0 then (startPoint = newPoint)
          ) else if (candidateValue > maxValues_i) then (
            maxValues = replace(i, candidateValue, maxValues)
          ); 
        );
      )
    );

  -- Take ceiling for all values, at every coordinate
  startPoint = apply(startPoint, value -> ceiling(value));
  minValues = apply(minValues, value -> ceiling(value));
  maxValues = apply(maxValues, value -> ceiling(value));
  
  return {startPoint, minValues, maxValues}
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that creates strings with varying lengths of tab indents

-- Appears in most print-related lines of code throughout the program

-- input: count, a positive int representing the number of tab indents desired

-- output: a string that has count number of tab spaces
--------------------------------------------------------------------------------

tabMult = (count) -> 
(
  theTabs := "";

  for i from 1 to count do(
    theTabs = concatenate(theTabs, "\t")
  );

  return theTabs
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that computes halfspace inequalties which describes the Newton 
--   Polyhedron for the rational power of an ideal.

-- Appears in rPowerStaircase function

-- input: an ideal I and a rational number r

-- output: a sequence with two elements. Each element is a matrix. The first 
--           matrix represents the coefficients of the variables, the second 
--           matrix represents the constants of the inequalities.
--------------------------------------------------------------------------------

getHalfSpaces = (I, r) ->
(
  -- Get NP(I)
  N := Npolyhedron I; 

  -- Get halfspace inequalities
  hspaces := halfspaces N; 

  -- Scale halfspaces by multiplying the constants only (i.e. second matrix)
  hspaces = (hspaces_0, r * hspaces_1); 

  -- Remove "trivial" inequalities, x_i >=0. Saves on caclculations
  zeroConstants := {};

  for i from 0 to numRows(hspaces_1) - 1 do(
    if hspaces_1_0_i == 0 then (
      zeroConstants = append(zeroConstants, i)
    );
  );

  hspaces = (submatrix'(hspaces_0, zeroConstants,), submatrix'(hspaces_1, zeroConstants,));

  return hspaces
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that checks if a point is in or outside of the Newton Polyhedron

-- Appears in actionHere function

-- input: a list (currPoint) representing a point and a sequence (hspaces)
--          representing the halfspace inequalities.

-- output: a boolean, true if the point is in the NP, false if it is not in NP
--------------------------------------------------------------------------------

checkInPolyhedron = (currPoint, hspaces) ->
(
  -- transpose currPoint so that we can multiply it to hspaces in next step
  matrixCoords := transpose matrix{currPoint}; 

  -- multiplication of point to coefficients of hspace inequaltiies
  result := hspaces_0 * matrixCoords; 

  -- some 'tracking' variables for while loop below
  rowCount := numRows result;
  n := 0;
  inPolyhedron := true;

  -- loop checks if inequalties are all satisfied
  while n < rowCount and inPolyhedron do (
    inPolyhedron = inPolyhedron and result_0_n <= (hspaces_1)_0_n;
    n = n + 1;
  );

  return inPolyhedron
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that moves a point by one unit in the relevant direction

-- Appears in actionHere fucntion

-- input: point, a list representing a point, 
--        direction, a boolean representing the direction of movement of the 
--          point, and 
--        varIndex, a positive int (varIndex) representing the coordinate which 
--          we are moving on. true  goes in the positive direction, while false 
--          goes in the negative direction.

-- output: a list representing the point that results from the move
--------------------------------------------------------------------------------

movePoint = (point, direction, varIndex) -> 
(
  i := varIndex - 1;
  newValue := if direction then point_i + 1 else point_i - 1;

  -- constructing a new list that represents the point resulting from the move
  newPoint := replace(i, newValue, point);

  return newPoint
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that checks if a point is within the boundary of the staircase 
--   algorithm. Points on the boundary itself is considered to be within.

-- Appears in actionHere and getMinGens functions

-- input: a list (currPoint) representing a point and a list of lists 
--          (extremeValues) which contains information on the boundary values 
--          of the algorithm. The first and second elements of extremeValues  
--          are lists of integers that store the min and max values at each 
--          coordinate.

-- output: a boolean true if input currPoint is within the boundaries, false 
--           otherwise
--------------------------------------------------------------------------------

withinExtremes = (currPoint, extremeValues) -> 
(
  toReturn := true;
  i := 0;

  -- loop over every coordinate and check if within boundaries
  while toReturn and i < numgens(R) do (
    toReturn = currPoint_i >= extremeValues_0_i 
                 and currPoint_i <= extremeValues_1_i;
    i = i + 1;
  );

  return toReturn
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that checks if a point is exactly on boundary of its traversal. 
--   The boundary value is contextual, it depends on the direction of 
--   traversal. It can either the min or the max of the input coordinate + 1.
--
--   To determine this, the function uses an extra argument isIncreasing
--   that tells it whether the point is trying to minimize or maximize the 
--   value at the currVarDim coordinate in the algorithm's traversal. If it  
--   seeks to minimize the input coordinate, the boundary value is the max of 
--   the input coordinate + 1. If it seeks to maximize the input coordinate, 
--   the boundary value is the min of the input coordinate + 1.

-- Appears in actionHere

-- input: currPoint, a list representing a point, 
--        currVarDim, a positive int representing the coordinate of interest, 
--        extremeValues, a list of lists which contains information on the 
--          boundary values of the algorithm,
--        isIncreasing, a boolean which tells us if it is true that the 
--          coordinate of interest is being maximized

-- output: a boolean that tells us whether the point is exactly on the boundary 
--           of its traversal, based on the explanation earlier.
--------------------------------------------------------------------------------

checkOnExtreme = (currPoint, currVarDim, extremeValues, isIncreasing) -> 
(
  -- minMax variable basically tells us which of the min or max we want, based 
  --   on the explanation above.

  minMax := if isIncreasing then 0 else 1;

  return {currPoint_currVarDim == extremeValues_minMax_currVarDim, minMax}
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that moves a point that is inside the NP to some new point near the boundary. This new point will only be diff in the var + 1 coordinate and will serve as the starting point to the waterfall that comes after it.

-- Appears in waterfall function

-- input: currPoint, a list representing a point, 
--        currVarDim, a positive int representing the coordinate of interest, 
--        extremeValues, a list of lists which contains information on the 
--          boundary values of the algorithm,
--        isIncreasing, a boolean which tells us if it is true that the 
--          coordinate of interest is being maximized
--        hspaces, a sequence representing the halfspace inequalities, 

-- output: a list representing a point that "just crosses" the boundary of the NP, by moving only parallel to the next coordinate axis.
--------------------------------------------------------------------------------

moveStartPointToBoundary = (currPoint, currInNP, currVarDim, extremeValues, hspaces) -> 
(
  -- direction in which we move the point is in the next variable direction, which depends on isIncreasing
  nextVarDirection := not currInNP;
  
  -- record which side of the NP we began this process, we want to "just cross" the boundary, into the opposite side from where we started
  origSide := currInNP;

  -- bounds to ensure that the new start point won't be beyond our extremes
  lowerBound := extremeValues_0_currVarDim;
  upperBound := extremeValues_1_currVarDim;

  while origSide == currInNP do (
    nextPoint := movePoint(currPoint, nextVarDirection, currVarDim + 1);

    -- if nextPoint falls outside the boundaries, we don't want it
    -- the if condition can certainly be optimized
    if lowerBound <= nextPoint_currVarDim and nextPoint_currVarDim <=  upperBound then (
      currPoint = nextPoint;
      currInNP = checkInPolyhedron(nextPoint, hspaces);
    ) else (
      break
    );
  );

  return {currPoint, currInNP}
);


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that implements the waterfall step (dimensional reduction) of the 
--   algorithm.

-- Appears in actionHere function

-- input: currPoint, a list representing a point, 
--        currVarDim, a positive int representing the coordinate of interest, 
--        extremeValues, a list of lists which contains information on the 
--          boundary values of the algorithm,
--        isIncreasing, a boolean which tells us if it is true that the 
--          coordinate of interest is being maximized
--        hspaces, a sequence representing the halfspace inequalities, 

-- output: a list of lists, representing the generators obtained in the 
--           waterfall step
--------------------------------------------------------------------------------

waterfall = (currPoint, currInNP, currVarDim, extremeValues, hspaces, isIncreasing) -> 
(
  -- initialize
  waterfallGens := {};

  -- waterfall is performed only at 3-dimensions and above, which is to say, 
  --   there must be at least two reminaing dimensions for a waterfall to make 
  --   sense. If 2-D, a belt step alone will suffice, that is accomplished in 
  --   actionHere. 

  if currVarDim <= numgens(R) - 2  then (
    print concatenate(tabMult(currVarDim - 1), "-- waterfall start");

    -- (hopefully) ensures that startpoint of this traversal is some point near the edge of NP
    if currInNP then (
      left := movePoint(currPoint, false, currVarDim + 1);
      bottom := movePoint(currPoint, false, currVarDim + 2);

      if checkInPolyhedron(left, hspaces) and checkInPolyhedron(bottom, hspaces) then (
        resultIn := moveStartPointToBoundary(currPoint, currInNP, currVarDim + 1, extremeValues, hspaces);
        currPoint = resultIn_0;
        currInNP = resultIn_1;
      )
    ) else (
      up := movePoint(currPoint, true, currVarDim + 2);
      right := movePoint(currPoint, true, currVarDim + 1);

      if checkInPolyhedron(up, hspaces) and checkInPolyhedron(right, hspaces) then (
        resultOut := moveStartPointToBoundary(currPoint, currInNP, currVarDim + 1, extremeValues, hspaces);
        currPoint = resultOut_0;
        currInNP = resultOut_1;
      )
    );

    -- currPoint = moveStartPointToBoundary(currPoint, currInNP, currVarDim + 1, extremeValues, hspaces);
    print concatenate(tabMult(currVarDim - 1), "currPoint moved to ", toString(currPoint));


    -- In 3-D, below would be the "moving down" case, where z-coordinate is decreasing but y-coordinate is increasing.
    print concatenate(tabMult(currVarDim - 1), "   increasing Dim #", toString(currVarDim + 1));
    minGensFromIncreasing := getMinGens(currPoint, currInNP, currVarDim + 1, extremeValues, hspaces, true);

    -- In 3-D, below would be the "moving up" case, where z-coordinate is increasing but y-coordinate is decreasing.
    print concatenate(tabMult(currVarDim - 1), "   decreasing Dim #", toString(currVarDim + 1));
    minGensFromDecreasing := getMinGens(currPoint, currInNP, currVarDim + 1, extremeValues, hspaces, false);

    -- join the lists of generators obtained from the two waterfall directions 
    --   and combine them into one list.
    waterfallGens = minGensFromIncreasing | minGensFromDecreasing;
    waterfallGens = unique(waterfallGens);

    if length(waterfallGens) != 0 then (
      print concatenate(, tabMult(currVarDim - 1), "-- waterfall end, generators collected: ", toString(waterfallGens));
    );
  );

  return waterfallGens
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that performs the belt step of the algorithm

-- Appears in getMinGens function

-- input: currPoint, a list representing a point, 
--        currVarDim, a positive int representing the coordinate of interest, 
--        extremeValues, a list of lists which contains information on the 
--          boundary values of the algorithm,
--        isIncreasing, a boolean which tells us if it is true that the 
--          coordinate of interest is being maximized
--        hspaces, a sequence representing the halfspace inequalities, 
--        previouslyIn, a boolean which tells us whether the previous point of 
--          the algorithm was within the NP.

-- output: a list of four elements.
--         The first is a list representing the next point that is to be 
--           visited in the algorithm,
--         the second is a list of lists representing the generators, through 
--           the algorithm, by visiting this point on the belt,
--         the third is a boolean that records whether the current point was 
--           inside NP,
--         the fourth is a boolean telling us if the next point is in NP.
--------------------------------------------------------------------------------

actionHere = (currPoint, currInNP, currVarDim, extremeValues, hspaces, isIncreasing, previouslyIn) -> 
( 
  nextPoint := {};
  newGens := {};

  extremCheck := checkOnExtreme(currPoint, currVarDim, extremeValues, isIncreasing);
  isOnExtreme := extremCheck_0;

  location := if currInNP then " in NP" else " outside NP";
  poleStr := if extremCheck_1 == 0 then "min" else "max";
  onExtremeStr := if isOnExtreme then concatenate(" and on ", poleStr, " extreme at Dim #", toString(currVarDim + 1)) else "";  
  print concatenate(tabMult(currVarDim - 1), toString(currPoint), location, onExtremeStr);


  -- direction for when we are outside NP
  currVarDirection := isIncreasing;

  -- direction for when we are inside NP
  nextVarDirection := not isIncreasing; 

  -- If the move prior to this was in the direction of the current variable, we 
  --   have to perform a waterfall step. To determine if the move was in the 
  --   current variable direction, it suffices to check if previouslyIn is the 
  --   opposite from isIncreasing.
  -- Moreover, if we are already on the extreme (note: we're still within NP, 
  --   as determined by getMinGens function called before entering actionHere), 
  --   then we need to perform a waterfall because the only moves we should  
  --   make are to move along the current variable direction and perpetually 
  --   get waterfalls to catch all other generators at the boundary, especially 
  --   those hiding in lower dimensions.
  if (previouslyIn != isIncreasing) or isOnExtreme then (
    newGens = newGens | waterfall(currPoint, currInNP, currVarDim, extremeValues, hspaces, isIncreasing);
  );

  if (isOnExtreme and (currVarDim <= numgens(R) - 2)) then (
      nextPoint = movePoint(currPoint, currVarDirection, currVarDim)
  ) else if not currInNP then (
    if isIncreasing then (
      nextPoint = movePoint(currPoint, currVarDirection, currVarDim)
    ) else (
      nextPoint = movePoint(currPoint, nextVarDirection, currVarDim + 1)
    ) 
  )else if currInNP then (
    if isIncreasing then (
      nextPoint = movePoint(currPoint, nextVarDirection, currVarDim + 1)
    ) else (
      nextPoint = movePoint(currPoint, currVarDirection, currVarDim)
    )
  );
  
  nextInNP := checkInPolyhedron(nextPoint, hspaces);
  nextPointOutside := not nextInNP;
  
  -- todo: check if this ocndition is truly necessary
  nextPointInsideButOB := nextInNP and not withinExtremes(nextPoint, extremeValues);
  
  
  -- adds current point as a generator if following conditions hold
  if currInNP and (nextPointOutside or nextPointInsideButOB) then (
    newGens = append(newGens, currPoint);
    newGens = unique(newGens);

    if nextPointOutside then (
      print concatenate(tabMult(currVarDim - 1), "generator added from exiting NP: ", toString(currPoint));
    ) else (
      print concatenate(tabMult(currVarDim - 1), "generator added from going OB: ", toString(currPoint));
    );
  );

  print concatenate(tabMult(currVarDim - 1), "done at this point\n");

  return {nextPoint, newGens, currInNP, nextInNP}
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that coordinates staircase algorithm, main helper called in this 
--   function is actionHere

-- Appears in rPowerStaircase function

-- input: currPoint, a list representing a point, 
--        currInNP, a boolean which tells is if currPoint is in NP
--        currVarDim, a positive int representing the coordinate of interest, 
--        extremeValues, a list of lists which contains information on the 
--          boundary values of the algorithm,
--        hspaces, a sequence representing the halfspace inequalities, 
--        isIncreasing, a boolean which tells us if it is true that the 
--          coordinate of interest is being maximized

-- output: a list of lists, where each list element represents a generator 
--           found through the staircase algorithm
--------------------------------------------------------------------------------

getMinGens = (currPoint, currInNP, currVarDim, extremeValues, hspaces, isIncreasing) -> 
(
  print concatenate(tabMult(currVarDim - 1), "Dim #", toString(currVarDim));

  -- initialize
  minGensArray := {}; 

  -- TODO:: call the 2-d algorithm for conciseness
  if currVarDim == numgens(R)  then (
    return minGensArray
  );

  -- variable that tells us if previous point viisted was in NP, initialized to 
  --   false so that a waterfall step can be initiated at the starting point
  previouslyIn := false; 

  -- while in NP continue performing algorithm
  while withinExtremes(currPoint, extremeValues) do (
    -- calls workhorse helper actionHere
     results := actionHere(currPoint, currInNP, currVarDim, extremeValues, hspaces, isIncreasing, previouslyIn);

    -- redefines variables for next loop
    currPoint = results_0;
    minGensArray = unique(minGensArray | results_1);
    previouslyIn = results_2;
    currInNP = results_3;
  ); 

  print concatenate(tabMult(currVarDim - 1), toString(currPoint), " out of bounds, ", "Dim #", toString(currVarDim), " traversal complete");
  print concatenate(tabMult(currVarDim - 1), "  * generators found at the end of belt traversal: ", toString(minGensArray), "\n");

  return minGensArray
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that makes an ideal out of list of points. These points represent 
--   the generators

-- Appears in rPowerStaircase

-- input: gensList, a list of lists where each list element represents a 
--          generator

-- output: a monomialIdeal constructed from the input list of generators
--------------------------------------------------------------------------------

makeIdeal = (gensList) ->
(
  -- INITIALIZE 
  monomialGens := {};

  -- LOOP THROUGH EACH GENERATOR AND VARIABLE
  for genIndex from 0 to (length(gensList) - 1) do (
    oneGenerator := {};

    for varIndex from 0 to (numgens(R) - 1) do (
      oneGenerator = append(oneGenerator, R_varIndex^((gensList_genIndex)_varIndex))
    );

    -- MULTIPLY THE VARIABLES USING FOLD, THEN ADD TO monomialGens array
    oneGenerator = fold((i, j) -> i * j, oneGenerator);
    monomialGens = append(monomialGens, oneGenerator);
  );

  -- CREATE THE IDEAL AND RETURN
  idealToReturn := monomialIdeal(monomialGens);

  return idealToReturn
);



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Function that is to be called when making the rational power of an ideal 

-- input: a monomial ideal I and a rational number r

-- output: a monomialIdeal representing rational power of I
--------------------------------------------------------------------------------

rPowerStaircase = (I,r) -> 
(
  -- get starting point and boundary values
  extremes := getExtremeValues(I, r);

  currPoint := extremes_0;
  extremeValues := drop(extremes, 1);

  -- get halfspaces
  hspaces := getHalfSpaces(I, r);

  -- get generators through staircase algorithm
  currInNP := checkInPolyhedron(currPoint, hspaces);
  minGenList := getMinGens(currPoint, currInNP, 1, extremeValues, hspaces, true);

  -- make the ideal using generators found
  newIdeal := makeIdeal(minGenList);
  
  print newIdeal;
  print toString(newIdeal);

  print "\n---------staircase algorithm complete---------\n";
  
  return newIdeal
);



----------------------------------------------------------------------------------------------------------------------------------------------------------------

