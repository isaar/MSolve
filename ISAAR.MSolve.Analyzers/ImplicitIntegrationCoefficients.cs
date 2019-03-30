using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class ImplicitIntegrationCoefficients
    {
        public double Mass { get; set; } = double.NaN;
        public double Damping { get; set; } = double.NaN;
        public double Stiffness { get; set; } = double.NaN;
    }
}
