using ISSAR.MSolve.IGAPreProcessor;
using ISSAR.MSolve.IGAPreProcessor.Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGASamplesConsole
{
    class Program
    {
        private static void SolveIsogeometricModelExample()
        {
            IGAModel model = new IGAModel();
            IGAModelReader.CreateModelFromFile(model, "Cantilever2D biquadratic");
        }

        static void Main(string[] args)
        {
            SolveIsogeometricModelExample();
        }


    }
}
