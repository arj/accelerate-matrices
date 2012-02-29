-- | This module provides the definition of a sparse matrix
-- and a multiplication.
module Data.Array.Accelerate.Math.SMVM where

import Data.Array.Accelerate.Math.SMVM.Matrix
import Data.Array.Accelerate.Types

import Data.Array.Unboxed
import Data.Array.Accelerate           (Vector, Segments, Acc, Array)
import qualified Data.Array.Accelerate as Acc
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

-- | Sparse-matrix vector multiplication
smvmAcc :: AccSparseMatrix Float -> AccVector Float -> AccVector Float
smvmAcc (segd, (inds, vals)) vec
  = let
      vecVals  = Acc.backpermute (Acc.shape inds) (\i -> Acc.index1 $ inds Acc.! i) vec
      products = Acc.zipWith (*) vecVals vals
    in
    Acc.foldSeg (+) 0 products segd

-- | Sparse-matrix vector multiplication. Wrapper for non-acc input.
smvm2Acc :: SparseMatrix Float -> Vector Float -> Acc (Vector Float)
smvm2Acc (segd', (inds', vals')) vec'
  = let
      segd     = Acc.use segd'
      inds     = Acc.use inds'
      vals     = Acc.use vals'
      vec      = Acc.use vec'
    in
      smvmAcc (segd, (inds, vals)) vec

-- | Returns the row number of a given sparse matrix
-- It is impossible to know the number of columns
smrows :: Num a => SparseMatrix a -> Int
smrows (s, _) = Acc.arraySize (Acc.arrayShape s)

-- | Creates a unity matrix of size n*n
smunity :: Int -> SparseMatrix Float
smunity n = (segments, (vectors, values))
 where
   segments     = Acc.fromList (Acc.Z Acc.:. n) $ take n $ repeat 1
   vectors      = Acc.fromList (Acc.Z Acc.:. n) [1..n]
   values       = Acc.fromList (Acc.Z Acc.:. n) $ take n $ repeat 1.0
