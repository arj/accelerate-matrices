-- | This module provides some basic vector operations.
module Data.Array.Accelerate.Math  where

import Data.Array.Accelerate.Types
import Data.Array.Unboxed
import Data.Array.Accelerate as Acc

-- | Returns the inverse of the given vector
vectorInverseAcc :: (Elt a, IsFloating a) => AccVector a -> AccVector a
vectorInverseAcc v = Acc.map (1/) v

(*^) :: (Elt a, IsFloating a) => AccVector a -> AccVector a -> AccVector a
(*^) = Acc.zipWith (*)

(*.) :: (Elt a, IsFloating a) => AccScalar a -> AccVector a -> AccVector a
alpha *. v = Acc.map (the alpha*) v

(./) :: (Elt a, IsFloating a) => AccVector a -> AccScalar a -> AccVector a
v ./ alpha = Acc.map (the alpha/) v

(/.) :: (Elt a, IsFloating a, Shape sh) => AccScalar a -> Acc (Acc.Array sh a) -> Acc (Acc.Array sh a)
alpha /. v = Acc.map (the alpha/) v

-- | Dot product for two vectors. Wrapper for non-acc input.
dotp2Acc :: (Elt a, IsFloating a) => Vector a -> Vector a -> Acc (Scalar a)
dotp2Acc xs ys =
  let xs' = use xs
      ys' = use ys
  in
      dotpAcc xs' ys'

-- | Dot product for two vectors.
dotpAcc :: (Elt a, IsFloating a) => AccVector a -> AccVector a -> AccScalar a
dotpAcc xs ys = Acc.fold (+) 0 (Acc.zipWith (*) xs ys)
