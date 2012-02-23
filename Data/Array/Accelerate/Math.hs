-- | This module provides some basic vector operations.
module Data.Array.Accelerate.Math  where

import Data.Array.Accelerate.Types
import Data.Array.Unboxed
import Data.Array.Accelerate as Acc

-- | Returns the inverse of the given vector
vectorInvertAcc :: AccVector Float -> AccVector Float
vectorInvertAcc v = Acc.map (1/) v

-- | Dot product for two vectors. Wrapper for non-acc input.
dotp2Acc :: Vector Float -> Vector Float -> Acc (Scalar Float)
dotp2Acc xs ys
  = let
      xs' = use xs
      ys' = use ys
    in
      dotpAcc xs' ys'

-- | Dot product for two vectors.
dotpAcc :: AccVector Float -> AccVector Float -> AccScalar Float
dotpAcc xs ys = Acc.fold (+) 0 (Acc.zipWith (*) xs ys)
