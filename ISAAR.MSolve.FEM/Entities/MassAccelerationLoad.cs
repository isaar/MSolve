using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.FEM.Entities
{
    public class MassAccelerationLoad
    {
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}
