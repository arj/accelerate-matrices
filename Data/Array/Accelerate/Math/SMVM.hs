{-#LANGUAGE ExplicitForAll,ScopedTypeVariables #-}
-- | This module provides the definition of a sparse matrix
-- and a multiplication.
module Data.Array.Accelerate.Math.SMVM where

import Data.Array.Accelerate.Math.SMVM.Matrix
import Data.Array.Accelerate.Types

import Data.Array.Unboxed as AU
--import Data.Array.Accelerate           (Vector, Segments, Acc, Array, use, Z, (:.))
import Data.Array.Accelerate
import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Array.Sugar as Sugar
import qualified Data.Vector.Unboxed   as V

-- | Definition of a sparse vector, i.e. a tuple with the first
-- entry beeing the column of the entry and the second beeing the
-- entry itself.
type SparseVector a = (Vector Int, Vector a)

-- | Definition of a sparse matrix based on the condensed-row format, i.e.
-- Segments as the number of non-zero element per row, and a 'SparseVector'
-- covering the non-zero entries.
type SparseMatrix a = (Segments Int, SparseVector a)

-- | Wrapper for 'SparseVector' type.
type AccSparseVector a = ((AccVector Int), (AccVector a))

-- | Wrapper for 'SparseMatrix' type.
type AccSparseMatrix a = ((AccSegments Int), (AccSparseVector a))

-- | Transfers a sparse matrix to accelerate by running use for all components.
usesm :: (Elt a, IsFloating a) => SparseMatrix a -> AccSparseMatrix a
usesm (segments, (vectors, values)) = (use segments, (use vectors, use values))

-- | Sparse-matrix vector multiplication
smvmAcc :: AccSparseMatrix Float -> AccVector Float -> AccVector Float
smvmAcc (segd, (inds, vals)) vec
  = let
      vecVals  = backpermute (shape inds) (\i -> index1 $ inds A.! i) vec
      products = A.zipWith (*) vecVals vals
    in
    foldSeg (+) 0 products segd

-- | Sparse-matrix vector multiplication. Wrapper for non-acc input.
smvm2Acc :: SparseMatrix Float -> Vector Float -> Acc (Vector Float)
smvm2Acc (segd', (inds', vals')) vec'
  = let
      segd     = use segd'
      inds     = use inds'
      vals     = use vals'
      vec      = use vec'
    in
      smvmAcc (segd, (inds, vals)) vec

-- | Builds a sparse matrix from an Accelerate 'A.Array' with an arbitrary sparse element.
fromArray :: (Eq a, Elt a, IsFloating a) => a -> A.Array DIM2 a -> SparseMatrix a
fromArray zero arr = toArraySM $ transform arr
  where
    empty   = ([], ([], []))
    (_ :. row_min :. col_min, _ :. row_max :. col_max) = Sugar.shapeToRange $ arrayShape arr
    rows_i  = [row_min..row_max]
    cols_i  = [col_min..col_max]
    toArray l = fromList (Z :. length l) l
    toArraySM (segList, (vecList, valList)) = (toArray segList, (toArray vecList, toArray valList))
    --
    innerTransform arr r = foldr (\c (vec, val) ->
                                            let e = indexArray arr (Z :. r :. c) in
                                            if e == zero then
                                              (vec, val)
                                            else
                                              (c : vec, e : val)

                                          ) ([],[]) cols_i

    --
    transform arr = foldr (\r (seg, (vec, val)) ->
                            let (vec1, val1) = innerTransform arr r
                            in
                             (length vec1 : seg, (vec1 ++ vec, val1 ++ val))

                      ) empty rows_i

-- | Builds a sparse matrix from an Accelerate 'A.Array' with 0.0 as the
-- sparse element.
fromArrayZero :: forall a.(Eq a, Elt a, IsFloating a, Fractional a) => A.Array DIM2 a -> SparseMatrix a
fromArrayZero = fromArray (0.0 :: a)


-- | Returns the row number of a given sparse matrix
-- It is impossible to know the number of columns
smrows :: Num a => SparseMatrix a -> Int
smrows (s, _) = arraySize (arrayShape s)

-- | Calculates the row number of an accelerate wrapped sparse matrix.
-- It is impossible to know the number of columns.
smrowsAcc :: Num a => AccSparseMatrix a -> Exp Int
smrowsAcc (s, _) = size s

-- TODO
-- smmul :: (Num a) => SparseMatrix a -> SparseMatrix a -> SparseVector a
-- smmul (seg1, (vec1, val1), (row1, col1)) (seg2, (vec2, val2), (row2, col2)) =
  


-- | Creates a unity matrix of size n*n
smunity :: Int -> SparseMatrix Float
smunity n = (segments, (vectors, values))
 where
   segments     = fromList (Z :. n) $ take n $ repeat 1
   vectors      = fromList (Z :. n) [1..n]
   values       = fromList (Z :. n) $ take n $ repeat 1.0

-- | Creates a unity matrix to be used in Accelerate.
smunity2Acc :: Int -> AccSparseMatrix Float
smunity2Acc n = (segments, (vectors, values))
 where
   segments     = use $ fromList (Z :. n) $ take n $ repeat 1
   vectors      = use $ fromList (Z :. n) [1..n]
   values       = use $ fromList (Z :. n) $ take n $ repeat 1.0
