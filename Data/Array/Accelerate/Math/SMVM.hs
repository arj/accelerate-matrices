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
import Prelude hiding (replicate, zip, unzip, map, scanl, scanl1, scanr, scanr1, zipWith,
                         filter, max, min, not, fst, snd, curry, uncurry, sum, head, tail,
                         drop, take, null, length, reverse, init, last, product, minimum,
                         maximum)
import qualified Prelude

-- | Definition of a sparse vector, i.e. a tuple with the first
-- entry beeing the column of the entry, the second beeing the
-- entry itself, and the third entry is the original size of the
-- vector.
--
type SparseVector a = (Vector Int, Vector a, Int)

-- | Definition of a sparse matrix based on the condensed-row format, i.e.
-- Segments as the number of non-zero element per row, and a 'SparseVector'
-- (without the size) covering the non-zero entries. The last entry
-- is the number of columns from the original matrix. The number of
-- rows can be extracted by length segments.
--
type SparseMatrix a = (Segments, (Vector Int, Vector a), Int)

-- | Wrapper for 'SparseVector' type.
--
type AccSparseVector a = (AccVector Int, AccVector a, AccScalar Int)

-- | Wrapper for 'SparseMatrix' type.
--
type AccSparseMatrix a = (AccSegments, (AccVector Int, AccVector a), AccScalar Int)

-- | Transfers a sparse matrix to accelerate by running use for all components.
usesm :: (Elt a, IsFloating a) => SparseMatrix a -> AccSparseMatrix a
usesm (segments, (vectors, values), size) = (use segments, (use vectors, use values), unit (constant size))

-- | Sparse-matrix vector multiplication
smvmAcc :: AccSparseMatrix Float -> AccVector Float -> AccVector Float
smvmAcc (segd, (inds, vals), _) vec
  = let
      vecVals  = backpermute (shape inds) (\i -> index1 $ inds A.! i) vec
      products = A.zipWith (*) vecVals vals
    in
    foldSeg (+) 0 products segd

-- | Sparse-matrix vector multiplication. Wrapper for non-acc input.
smvm2Acc :: SparseMatrix Float -> Vector Float -> Acc (Vector Float)
smvm2Acc (segd', (inds', vals'), cols') vec'
  = let
      segd     = use segd'
      inds     = use inds'
      vals     = use vals'
      vec      = use vec'
      cols     = unit $ constant cols'
    in
      smvmAcc (segd, (inds, vals), cols) vec

-- | Builds a sparse matrix from an Accelerate 'A.Array' with an arbitrary sparse element.
fromArray :: (Eq a, Elt a, IsFloating a) => a -> A.Array DIM2 a -> SparseMatrix a
fromArray zero arr = toArraySM $ transform arr
  where
    empty   = ([], ([], []),0)
    (_ :. row_min :. col_min, _ :. row_max :. col_max) = Sugar.shapeToRange $ arrayShape arr
    rows_i  = [row_min..row_max]
    cols_i  = [col_min..col_max]
    toArray l = fromList (Z :. Prelude.length l) l
    toArraySM (segList, (vecList, valList), cols) = (toArray segList, (toArray vecList, toArray valList), cols)
    --
    innerTransform arr r = foldr (\c (vec, val) ->
                                            let e = indexArray arr (Z :. r :. c) in
                                            if e == zero then
                                              (vec, val)
                                            else
                                              (c : vec, e : val)

                                          ) ([],[]) cols_i

    --
    transform arr = foldr (\r (seg, (vec, val), _) ->
                            let (vec1, val1) = innerTransform arr r
                            in
                             (Prelude.length vec1 : seg, (vec1 ++ vec, val1 ++ val), (col_max - col_min + 1))

                      ) empty rows_i

-- | Builds a sparse matrix from an Accelerate 'A.Array' with 0.0 as the
-- sparse element.
fromArrayZero :: forall a.(Eq a, Elt a, IsFloating a, Fractional a) => A.Array DIM2 a -> SparseMatrix a
fromArrayZero = fromArray (0.0 :: a)

-- | Returns the row count of a given sparse matrix
--
smrows :: Num a => SparseMatrix a -> Int
smrows (s, _, _) = arraySize (arrayShape s)

-- | Calculates the row count of an accelerate wrapped sparse matrix.
--
smrowsAcc :: Num a => AccSparseMatrix a -> Exp Int
smrowsAcc (s, _, _) = size s

-- | Returns the col count of a given sparse matrix
--
smcols :: Num a => SparseMatrix a -> Int
smcols (_, _, cols) = cols

-- | Calculates the col count of an accelerate wrapped sparse matrix.
--
smcolsAcc :: Num a => AccSparseMatrix a -> Exp Int
smcolsAcc (_, _, cols) = the cols

-- TODO
-- smmul :: (Num a) => SparseMatrix a -> SparseMatrix a -> SparseVector a
-- smmul (seg1, (vec1, val1), (row1, col1)) (seg2, (vec2, val2), (row2, col2)) =
  


-- | Creates a unity matrix of size n*n
smunity :: Int -> SparseMatrix Float
smunity n = (segments, (vectors, values), n)
 where
   segments     = fromList (Z :. n) $ Prelude.take n $ repeat 1
   vectors      = fromList (Z :. n) [1..n]
   values       = fromList (Z :. n) $ Prelude.take n $ repeat 1.0

-- | Creates a unity matrix to be used in Accelerate.
smunity2Acc :: Int -> AccSparseMatrix Float
smunity2Acc n = (segments, (vectors, values), unit $ constant n)
 where
   segments     = use $ fromList (Z :. n) $ Prelude.take n $ repeat 1
   vectors      = use $ fromList (Z :. n) [1..n]
   values       = use $ fromList (Z :. n) $ Prelude.take n $ repeat 1.0


------------------------------


sparsedotpAcc :: AccSparseVector Float -> AccSparseVector Float -> AccScalar Float
sparsedotpAcc (idx1,val1,_) (idx2,val2,_) = A.foldAll (+) 0 $ A.map mapfun m1
  where
    v1 = A.zip idx1 val1 :: AccVector (Int, Float)
    v2 = A.zip idx2 val2 :: AccVector (Int, Float)
    genfun ix = let Z :. i :. j = unlift ix in lift (v1 A.! index1 i, v2 A.! index1 j)
    --
    m1 = generate (lift (Z :. size v1 :. size v2)) genfun
    --
    mapfun :: Exp ((Int, Float), (Int, Float)) -> Exp Float
    mapfun arg = let (ii,jj) = unlift arg :: (Exp (Int, Float), Exp (Int, Float)) in
                 let (i, vi) = unlift ii  :: (Exp Int, Exp Float) in
                 let (j, wj) = unlift jj  :: (Exp Int, Exp Float) in
                 (i ==* j) ? (vi * wj, constant 0)

