using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using ISAAR.MSolve.Geometry.Coordinates;

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

        public FixedPropagator(string anglesLengthsPath, Propagator actualPropagator = null)
        {
            this.actualPropagator = actualPropagator;
            checkPropagation = (actualPropagator == null) ? false : true;
            Logger = ReadFromFile(anglesLengthsPath);
            iteration = 0;
        }

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

        //public FixedPropagator CreateFromPath(IReadOnlyList<CartesianPoint> knownCrackPath, Propagator actualPropagator = null)
        //{
        //    int growthSteps = knownCrackPath.Count - 2; // The first 2 vertices are the initial crack.
        //    double[] angles = new double[growthSteps];
        //    double[] lengths = new double[growthSteps];
        //    for (int i = 0; i < growthSteps; ++i)
        //    {
        //        CartesianPoint previousTip = knownCrackPath[i];
        //        CartesianPoint currentTip = knownCrackPath[i+1];
        //        CartesianPoint nextTip = knownCrackPath[i+2];
        //        var oldSegment = new DirectedSegment2D(previousTip, currentTip);
        //        var newSegment = new DirectedSegment2D(currentTip, nextTip);
        //        lengths[i] = newSegment.Length;
        //        CartesianPoint local = oldSegment.TransformGlobalToLocalPoint(nextTip);
        //        angles[i] = Math.Atan2(local.Y, local.X);
        //    }
        //    return new FixedPropagator(angles, lengths, actualPropagator);
        //}

        
        public (double growthAngle, double growthLength) Propagate(Dictionary<int, Vector> totalFreeDisplacements, 
            CartesianPoint crackTip, TipCoordinateSystem tipSystem, IReadOnlyList<XContinuumElement2D> tipElements)
        {
            if (iteration >= Logger.GrowthLengths.Count) throw new IndexOutOfRangeException(
                $"Only {Logger.GrowthLengths.Count} iterations have been recorder.");
            double angle = Logger.GrowthAngles[iteration];
            double length = Logger.GrowthLengths[iteration];
            if (checkPropagation)
            {
                actualPropagator.Propagate(totalFreeDisplacements, crackTip, tipSystem, tipElements);
                Console.Write($"Growth angle: expected = {angle}");
                Console.WriteLine($"   -   computed = {actualPropagator.Logger.GrowthAngles[iteration]}");
                Console.Write($"Growth length: expected = {length}");
                Console.WriteLine($"   -   computed = {actualPropagator.Logger.GrowthLengths[iteration]}");
            }

            ++iteration;
            return (angle, length);
        }

        /// <summary>
        /// The file must have the number of iteration in the first line. Then in each line it must have the growth angle in the
        /// local cartesian system around the crack tip, the growth length SIF Mode I and SIF mode II of a new iteration, 
        /// separated by a space. No empty lines at the end.
        /// </summary>
        /// <param name="anglesLengthsPath"></param>
        /// <returns></returns>
        private PropagationLogger ReadFromFile(string anglesLengthsPath)
        {
            using (var reader = new StreamReader(anglesLengthsPath))
            {
                var logger = new PropagationLogger();
                int numIterations = int.Parse(reader.ReadLine());
                for (int i = 0; i < numIterations; ++i)
                {
                    string[] words = reader.ReadLine().Split(' ');
                    if (words.Length != 4) throw new IOException("Each line must have the growth angle, growth length," 
                        +" SIF Mode I and SIF Mode II of an iteration, separated by a space");
                    logger.GrowthAngles.Add(double.Parse(words[0]));
                    logger.GrowthLengths.Add(double.Parse(words[1]));
                    logger.SIFsMode1.Add(double.Parse(words[2]));
                    logger.SIFsMode2.Add(double.Parse(words[3]));
                }
                return logger;
            }
        }
    }
}
