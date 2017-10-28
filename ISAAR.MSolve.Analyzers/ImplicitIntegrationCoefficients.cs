using ISAAR.MSolve.Analyzers.Interfaces;

namespace ISAAR.MSolve.Analyzers
{
    public class ImplicitIntegrationCoefficients : IImplicitIntegrationCoefficients
    {
        public double Mass { get; set; }
        public double Damping { get; set; }
        public double Stiffness { get; set; }
    }
}
