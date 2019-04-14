using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.PreconditionedConjugateGradient
{
    /// <summary>
    /// Manages the insertion and removal of PCG direction vectors and related data, that will be used for reorthogonalization.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PcgReorthogonalizationCache
    {
        /// <summary>
        /// The conjugate direction vectors stored so far.
        /// </summary>
        public List<IVectorView> Directions { get; } = new List<IVectorView>();

        /// <summary>
        /// The products <see cref="Directions"/> * systemMatrix * <see cref="Directions"/> stored so far.
        /// </summary>
        public List<double> DirectionsTimesMatrixTimesDirections { get; } = new List<double>();

        /// <summary>
        /// The products systemMatrix * <see cref="Directions"/> stored so far.
        /// </summary>
        public List<IVectorView> MatrixTimesDirections { get; } = new List<IVectorView>();

        /// <summary>
        /// Discards the direction vectors and any corresponding data of the newest PCG iterations.
        /// </summary>
        /// <param name="numOldVectorsToRemove">
        /// The number of the newest entries (direction vectors and corresponding data) to discard. If it exceeds the number of
        /// entries currently stored, they will all be discarded without throwing any exceptions.
        /// </param>
        public void RemoveNewDirectionVectorData(int numNewVectorsToRemove)
        {
            if (numNewVectorsToRemove > Directions.Count)
            {
                Directions.Clear();
                MatrixTimesDirections.Clear();
                DirectionsTimesMatrixTimesDirections.Clear();
            }
            else
            {
                int start = Directions.Count - numNewVectorsToRemove;
                Directions.RemoveRange(start, numNewVectorsToRemove);
                MatrixTimesDirections.RemoveRange(start, numNewVectorsToRemove);
                DirectionsTimesMatrixTimesDirections.RemoveRange(start, numNewVectorsToRemove);
            }
        }

        /// <summary>
        /// Discards the direction vectors and any corresponding data of the oldest PCG iterations.
        /// </summary>
        /// <param name="numOldVectorsToRemove">
        /// The number of the oldest entries (direction vectors and corresponding data) to discard. If it exceeds the number of
        /// entries currently stored, they will all be discarded without throwing any exceptions.
        /// </param>
        public void RemoveOldDirectionVectorData(int numOldVectorsToRemove)
        {
            if (numOldVectorsToRemove > Directions.Count)
            {
                Directions.Clear();
                MatrixTimesDirections.Clear();
                DirectionsTimesMatrixTimesDirections.Clear();
            }
            else
            {
                Directions.RemoveRange(0, numOldVectorsToRemove);
                MatrixTimesDirections.RemoveRange(0, numOldVectorsToRemove);
                DirectionsTimesMatrixTimesDirections.RemoveRange(0, numOldVectorsToRemove);
            }
        }

        /// <summary>
        /// Stores a new direction vector and other related data. The new entries will be regarded as latest.
        /// </summary>
        /// <param name="pcg">The Preconditioned Conjugate Gradient Aglorithm that uses this object.</param>
        public void StoreDirectionData(PcgWithReorthogonalization pcg)
        {
            Directions.Add(pcg.Direction.Copy());
            MatrixTimesDirections.Add(pcg.MatrixTimesDirection.Copy());
            DirectionsTimesMatrixTimesDirections.Add(pcg.DirectionTimesMatrixTimesDirection);
        }
    }
}
