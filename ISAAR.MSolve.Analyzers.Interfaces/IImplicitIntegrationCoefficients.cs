using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IImplicitIntegrationCoefficients
    {
        double Mass { get; set; }
        double Damping { get; set; }
        double Stiffness { get; set; }
    }
}
