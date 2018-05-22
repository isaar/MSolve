using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.CrackPropagation
{
    /// <summary>
    /// Only for propagation from one tip for now.
    /// </summary>
    class FixedPropagator: IPropagator
    {
        /// <summary>
        /// In the local polar coordinate system defined at the crack tip.
        /// </summary>
        private readonly bool checkPropagation;
        private readonly Propagator actualPropagator;
        private int iteration;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="logger">Intermediate crack propagation data that were gathered by a previous analysis and will be 
        ///     enforced now.</param>
        /// <param name="actualPropagator">For asserting purposes</param>
        public FixedPropagator(PropagationLogger logger, Propagator actualPropagator = null)
        {
            this.actualPropagator = actualPropagator;
            checkPropagation = (actualPropagator == null) ? false : true;
            this.Logger = logger;
            for (int i = 0; i < logger.GrowthAngles.Count; ++i)
            {
                logger.InteractionIntegralsMode1.Add(0.0);
                logger.InteractionIntegralsMode2.Add(0.0);
                logger.SIFsMode1.Add(0.0);
                logger.SIFsMode2.Add(0.0);
            }

            iteration = 0;
        }

        public PropagationLogger Logger { get; }

        //public FixedPropagator CreateFromPath(IReadOnlyList<ICartesianPoint2D> knownCrackPath, Propagator actualPropagator = null)
        //{
        //    int growthSteps = knownCrackPath.Count - 2; // The first 2 vertices are the initial crack.
        //    double[] angles = new double[growthSteps];
        //    double[] lengths = new double[growthSteps];
        //    for (int i = 0; i < growthSteps; ++i)
        //    {
        //        ICartesianPoint2D previousTip = knownCrackPath[i];
        //        ICartesianPoint2D currentTip = knownCrackPath[i+1];
        //        ICartesianPoint2D nextTip = knownCrackPath[i+2];
        //        var oldSegment = new DirectedSegment2D(previousTip, currentTip);
        //        var newSegment = new DirectedSegment2D(currentTip, nextTip);
        //        lengths[i] = newSegment.Length;
        //        ICartesianPoint2D local = oldSegment.TransformGlobalToLocalPoint(nextTip);
        //        angles[i] = Math.Atan2(local.Y, local.X);
        //    }
        //    return new FixedPropagator(angles, lengths, actualPropagator);
        //}

        public (double growthAngle, double growthLength) Propagate(IDofOrderer dofOrderer, Vector totalFreeDisplacements,
            Vector totalConstrainedDisplacements)
        {
            double angle = Logger.GrowthAngles[iteration];
            double length = Logger.GrowthLengths[iteration];
            if (checkPropagation)
            {
                actualPropagator.Propagate(dofOrderer, totalFreeDisplacements, totalConstrainedDisplacements);
                Console.Write($"Growth angle: expected = {angle}");
                Console.WriteLine($"   -   computed = {actualPropagator.Logger.GrowthAngles[iteration]}");
                Console.Write($"Growth length: expected = {length}");
                Console.WriteLine($"   -   computed = {actualPropagator.Logger.GrowthLengths[iteration]}");
            }

            ++iteration;
            return (angle, length);
        }

    }
}
