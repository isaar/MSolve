using ISAAR.MSolve.Optimization.Commons;
using ISAAR.MSolve.Optimization.Problems;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Troschuetz.Random;

namespace ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public class StandardBinaryCoding : AbstractBinaryCoding
    {
        public StandardBinaryCoding(OptimizationProblem problem, int bitsPerContinuousVariable, int bitsPerIntegerVariable) :
                        base(problem, bitsPerContinuousVariable, bitsPerIntegerVariable)
        {
        }

        protected sealed override int BitstringToDecimalInteger(bool[] bits, int start, int length)
        {
            return BinaryUtilities.StandardBinaryToDecimal(bits, start, length);
        }

        protected override void DecimalIntegerToBitstring(int dec, bool[] bits, int start, int length)
        {
            BinaryUtilities.DecimalToStandardBinary(dec, bits, start, length);
        }
    }
}
