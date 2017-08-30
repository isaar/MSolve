using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation
{
    class PropagationLogger
    {
        public List<double> InteractionIntegralsMode1 { get; }
        public List<double> InteractionIntegralsMode2 { get; }
        public List<double> SIFsMode1 { get; }
        public List<double> SIFsMode2 { get; }
        public List<double> GrowthAngles { get; }
        public List<double> GrowthLengths { get; }

        public PropagationLogger()
        {
            InteractionIntegralsMode1 = new List<double>();
            InteractionIntegralsMode2 = new List<double>();
            SIFsMode1 = new List<double>();
            SIFsMode2 = new List<double>();
            GrowthAngles = new List<double>();
            GrowthLengths = new List<double>();
        }

        public void PrintAnalysisStep(int iteration)
        {
            Console.WriteLine("Analysis step " + iteration + ":");
            Console.WriteLine("Interaction integral - Mode 1 = " + InteractionIntegralsMode1[iteration]);
            Console.WriteLine("Interaction integral - Mode 2 = " + InteractionIntegralsMode2[iteration]);
            Console.WriteLine("SIF - Mode 1 = " + SIFsMode1[iteration]);
            Console.WriteLine("SIF - Mode 2 = " + SIFsMode2[iteration]);
            Console.WriteLine("Crack growth angle (local tip system) = " + GrowthAngles[iteration]);
            Console.WriteLine("Crack growth length = " + GrowthLengths[iteration]);
        }
    }
}
