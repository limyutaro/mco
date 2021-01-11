-- Modified staircase algorithm by iteratng over all k-2 variables, and concluding finally by a 2-d staircase traversal.
-- V1 JAN 7, 2020. 10:54

-- R = QQ[w,x,y,z];
-- I = ideal (x^7*y^3, x^5*z, y^5);
-- I = ideal(w*z, w*x, x^5);
-- r = 8/5;

-- rPowerStaircase(I,r)


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
loadPackage ("Polyhedra", Reload => true);

--------------------------------------------------------------------------------
getExtremeValuesStaircase = (I, r) ->
(
  generatorArray := apply(apply(flatten entries gens I, exponents), flatten);

  -- initialize minValues and maxValues to be the exponent values of the first generator
  minValues := r * generatorArray_0;
  maxValues := minValues;

  restOfGenerators := drop(generatorArray, 1);

  -- iterate through all generators to obtain the values that describe the hyperrectangle
  n := numgens(R)-1;
  
  scan(restOfGenerators, 
    point -> (
      newPoint := r * point;
      -- for each coordinate of the generator, check if it is the max or min at that coordinate
      for i from 0 to n do(
        candidateValue := newPoint_i;

        if (candidateValue < minValues_i) then (
          minValues = replace(i, candidateValue, minValues)
        ) else if (candidateValue > maxValues_i) then (
          maxValues = replace(i, candidateValue, maxValues)
        ); 
        
      );
    )
  );

  -- take ceilings to make them lattice points
  minValues = apply(minValues, value -> ceiling(value));
  maxValues = apply(maxValues, value -> ceiling(value));

  return {minValues, maxValues}
);

--------------------------------------------------------------------------------
NpolyhedronStaircase = I ->
(
  C := posHull id_(ZZ^(dim ring I));

  NP := C + newtonPolytope sum flatten entries gens I;

  return NP
);
--------------------------------------------------------------------------------
getHalfSpacesStaircase = (I, r) ->
(
  -- Get NP(I)
  N := NpolyhedronStaircase I; 

  -- Get halfspace inequalities
  hspaces := halfspaces N;

  -- Scale halfspaces by multiplying the constants only (i.e. second matrix)
  hspaces = (hspaces_0, r * hspaces_1); 

  -- Remove "trivial" inequalities, x_i >=0. Saves on caclculations
  zeroConstants := {};

  for i from 0 to numRows(hspaces_1) - 1 do(
    if hspaces_1_0_i == 0 then ( -- todo could potentially just change this to <0 condition
      zeroConstants = append(zeroConstants, i)
    );
  );

  hspaces = (submatrix'(hspaces_0, zeroConstants,), submatrix'(hspaces_1, zeroConstants,));

  return hspaces
);
--------------------------------------------------------------------------------
checkInPolyhedronStaircase = (point, hspaces) ->
(
  matrixCoords := transpose matrix{point}; -- convert point to a matrix x
  result := hspaces_0 * matrixCoords; -- compute Ax 

  numOfIneql := numRows result;
  i := 0;
  satisfiesIneql := true;

  -- CHECK IF INEQUALITIES ARE ALL SATISFIED
  while (i < numOfIneql) and satisfiesIneql do (
    satisfiesIneql = satisfiesIneql and result_0_i <= (hspaces_1)_0_i; --todo: can remove the satifiesIneql condition. We can just consider the latter logic check "result_0_i <= (hspaces_1)_0_i"
    i = i + 1;
  );

  return satisfiesIneql
);
--------------------------------------------------------------------------------
cascadeStaircase = (point, varNum, minVals, maxVals, hspaces) -> 
(
  generatorsList := {};
  n := numgens(R) - 1;

  if varNum < n - 1 then (
    for i from minVals_varNum to maxVals_varNum do(
      newpoint := replace(varNum, i, point);

      generatorsList = generatorsList | cascadeStaircase(newpoint, varNum + 1, minVals, maxVals, hspaces)
    )
  ) else if varNum == n - 1 then (
    generatorsList = generatorsList | staircaseStaircase(point, varNum, minVals, maxVals, hspaces)
  );

  return generatorsList
);
--------------------------------------------------------------------------------
staircaseStaircase = (currPoint, varNum, minVals, maxVals, hspaces) -> (
  --todo trim hspaces to essential inequalities here
  print("\nnew staircase");

  gensArray := {}; 

  x := varNum;
  y := varNum + 1;
  maxX := maxVals_x;
  minY := minVals_y;
  maxY := maxVals_y;
  
  -- boolean that tells us if the previously visited lattice point was in r * NP(I)
  previouslyIn := false; 

  withinHyperrect := true;

  -- performs 2-d staircase until we arrive at the maximum value for the last coordinate
  while withinHyperrect do (
    
    inPolyhedron := checkInPolyhedronStaircase(currPoint, hspaces);
    topLeft := currPoint_y == maxY and not inPolyhedron;

    if topLeft then ( -- point is outside but topleft
      print(concatenate(toString(currPoint), " in topleft, outside"));
      currPoint = replace(x, currPoint_x + 1, currPoint)
    ) else if inPolyhedron then ( -- point is inside

      bottomRight := currPoint_y == minY; 
      
      if bottomRight then ( -- check if on bottom right
        gensArray = append(gensArray, currPoint);
        print(concatenate(toString(currPoint), " in bottomright, inside"));
      ) else (
        print(concatenate(toString(currPoint), " inside"));
      );

      previouslyIn = true; 
      currPoint = replace(y, currPoint_y - 1, currPoint)
    
    ) else ( -- point is outside
      if previouslyIn then (
        gensArray = append(gensArray, replace(y, currPoint_y + 1, currPoint))
      );
      
      print(concatenate(toString(currPoint), " outside"));
      previouslyIn = false; 
      currPoint = replace(x, currPoint_x + 1, currPoint)
    );

    withinHyperrect = currPoint_x <= maxX and currPoint_y >= minY;
  );

  -- print("\n");
  return gensArray
)
--------------------------------------------------------------------------------
makeIdealStaircase = (gensArray) -> (
  monomialGens := {};

  -- loop through each generator and variable
  for genIndex from 0 to (length(gensArray) - 1) do (
    oneGenerator := {};

    for varIndex from 0 to (numgens(R) - 1) do (
      oneGenerator = append(oneGenerator, R_varIndex^((gensArray_genIndex)_varIndex))
    );

    -- multiply the variables using fold, then add to monomialGens array
    oneGenerator = fold((i, j) -> i * j, oneGenerator);
    monomialGens = append(monomialGens, oneGenerator);
  );

  idealToReturn := monomialIdeal(monomialGens);

  return idealToReturn
);
--------------------------------------------------------------------------------

rPowerStaircase = (I,r) -> (

  extremes := getExtremeValuesStaircase(I,r);
  minVals := extremes_0;
  maxVals := extremes_1;

  -- define point such that the (n-1)th coordinate is smallest and nth coordinate is largest. other coordinates are not important, because they will be changed by the cascade function anyway
  n := numgens(R) - 1;
  startPoint := replace(n, maxVals_n, minVals);

  hspaces := getHalfSpacesStaircase(I, r);

  minGens := cascadeStaircase(startPoint, 0, minVals, maxVals, hspaces);

  return makeIdealStaircase(minGens)
);

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


