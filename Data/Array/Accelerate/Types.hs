-- | This modules provides some basic shortcuts
--   for accelerate-constructor-wrapped types.
module Data.Array.Accelerate.Types where

import Data.Array.Accelerate

-- | Wrapper for 'Scalar' type.
type AccScalar a = Acc (Scalar a)

-- | Wrapper for 'Vector' type.
type AccVector a = Acc (Vector a)

-- | Wrapper for 'Segments' type.
type AccSegments a = Acc (Segments Int)
