using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Utilities;

namespace ISAAR.MSolve.XFEM.Tests.Geometry
{
    class AngleTests
    {
        public static void Main()
        {
            double perturbation = 1e-5;
            Console.WriteLine("[0, 2pi) \t -> \t (-pi, pi]");
            for (int i = 0; i < 12; ++i)
            {
                double angle = -3 * Math.PI + Math.PI / 2.0 * i;
                Console.WriteLine("{0} \t -> \t {1}", angle - perturbation, AngleUtilities.Wrap(angle - perturbation));
                Console.WriteLine("{0} \t -> \t {1}", angle, AngleUtilities.Wrap(angle));
                Console.WriteLine("{0} \t -> \t {1}", angle + perturbation, AngleUtilities.Wrap(angle + perturbation));
            }
        }
    }
}
