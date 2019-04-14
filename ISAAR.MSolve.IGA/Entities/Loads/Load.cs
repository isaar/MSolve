using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Load
    {
        public ControlPoint ControlPoint { get; set; }
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}
