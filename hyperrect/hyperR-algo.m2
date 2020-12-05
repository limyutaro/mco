-- Hyperrectangle RATIONAL POWER FINDER
-- program that will definitely find rational power. Only recommended for small ideals.
-- V2 AUG 9 14:30

-- changed name of algorihtm from naive to hyperrectangle, which corresponding changes to the name of the functions in this program



-- R = QQ[w,x,y,z];
-- I = ideal (x^7*y^3, x^5*z, y^5);
-- I = ideal(w*z, w*x, x^5);
-- r = 8/5;


-- rPowerHR(I,r)


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
loadPackage ("Polyhedra", Reload => true);

--------------------------------------------------------------------------------
getExtremeValuesHR = (I, r) ->
(
  generatorArray := apply(apply(flatten entries gens I, exponents), flatten);

  -- INITIALIZE LIST OF START AND END POINTS TO A LIST OF SIZE DIMENSIONS, EACH ELEMENT BEING THE FIRST GENERATOR IN generatorArray
  dummyGenerator := r * generatorArray_0;
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
        -- FOR EACH COORDINATE OF THE GENERATOR, CHECK IF IT IS THE MAX OR MIN AT THAT COORDINATE
        for i from 0 to (numgens(R) - 1) do(
          candidateValue := newPoint_i;

          if (candidateValue < minValues_i) then (
            minValues = replace(i, candidateValue, minValues)
          ) else if (candidateValue > maxValues_i) then (
            maxValues = replace(i, candidateValue, maxValues)
          ); 
        );
      )
    );

  -- TAKE A CEILING OF THEIR COORDINATES
  minValues = apply(minValues, value -> ceiling(value));
  maxValues = apply(maxValues, value -> ceiling(value));

  return {minValues, maxValues}
);

--------------------------------------------------------------------------------
NpolyhedronHR = I ->
(
  C := posHull id_(ZZ^(dim ring I));

  NP := C + newtonPolytope sum flatten entries gens I;

  return NP
);
--------------------------------------------------------------------------------
getHalfSpacesHR = (I, r) ->
(
  N := NpolyhedronHR I; -- GET NP(I)
  hspaces := halfspaces N; -- GET HALFSPACES
  hspaces = (hspaces_0, r * hspaces_1); -- SCALE HALFSPACES

  return hspaces
);
--------------------------------------------------------------------------------
checkInPolyhedronHR = (currPoint, hspaces) ->
(
  matrixCoords := transpose matrix{currPoint}; -- CONVERT currPoint TO MATRIX FORM
  result := hspaces_0 * matrixCoords; -- MULTIPLY POINT BY HALSPACE MATRIX

  -- SOME VARIABLES THAT HELP WITH WHILE LOOP AND TRACKER CALLED inPolyhedron
  rowCount := numRows result;
  n := 0;
  inPolyhedron := true;

  -- CHECK IF INEQUALITIES ARE ALL SATISFIED
  while n < rowCount and inPolyhedron do (
    inPolyhedron = inPolyhedron and result_0_n <= (hspaces_1)_0_n;
    n = n + 1;
  );

  return inPolyhedron
);
--------------------------------------------------------------------------------
makeIdealHR = (gensArray) ->
(
  -- INITIALIZE 
  monomialGens := {};

  -- LOOP THROUGH EACH GENERATOR AND VARIABLE
  for genIndex from 0 to (length(gensArray) - 1) do (
    oneGenerator := {};

    for varIndex from 0 to (numgens(R) - 1) do (
      oneGenerator = append(oneGenerator, R_varIndex^((gensArray_genIndex)_varIndex))
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


hyperRectHelper = (point, varNum, extremeVals, hspaces) -> 
(
  myList := {};

  if varNum < numgens(R) - 1 then (
    for i from (extremeVals_0)_varNum to (extremeVals_1)_varNum do(
      newpoint := replace(varNum, i, point);

      myList = myList | hyperRectHelper(newpoint, varNum + 1, extremeVals, hspaces)
    )
  ) else if varNum == numgens(R) - 1 then (
    for i from (extremeVals_0)_varNum to (extremeVals_1)_varNum do(
      newpoint := replace(varNum, i, point);
    
      if checkInPolyhedronHR(newpoint, hspaces) then(
        myList = append(myList, newpoint)
      )
    )
  );

  return myList
);

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

rPowerHR = (I,r) -> (
  extremes := getExtremeValuesHR(I,r);
  hspaces := getHalfSpacesHR(I, r);

  myList := hyperRectHelper(extremes_0, 0, extremes,hspaces);

  return makeIdealHR(myList)
);

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
