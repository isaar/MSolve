using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class ImplicitIntegrationCoefficients
    {
        public double Mass { get; set; }
        public double Damping { get; set; }
        public double Stiffness { get; set; }
    }
}
