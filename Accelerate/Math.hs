-- | This module provides some basic vector operations.
module Accelerate.Math  where

import Accelerate.Types
import Data.Array.Unboxed
import Data.Array.Accelerate as Acc

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
