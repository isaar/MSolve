using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.LinearAlgebra.Iterative.ConjugateGradient
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
        /// <param name="direction">
        /// The new direction vector, which is conjugate to other direction vectors stored previously.
        /// </param>
        /// <param name="matrixTimesDirection">
        /// The product systemMatrix * <paramref name="direction"/>.
        /// </param>
        /// <param name="directionTimesMatrixTimesDirection">
        /// The product <paramref name="direction"/> * systemMatrix * <paramref name="direction"/>.
        /// </param>
        public void StoreDirectionData(IVectorView direction, IVectorView matrixTimesDirection,
            double directionTimesMatrixTimesDirection)
        {
            Directions.Add(direction.Copy());
            MatrixTimesDirections.Add(matrixTimesDirection.Copy());
            DirectionsTimesMatrixTimesDirections.Add(directionTimesMatrixTimesDirection);
        }
    }
}
