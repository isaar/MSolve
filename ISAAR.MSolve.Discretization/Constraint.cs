using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.Discretization
{
    public class Constraint
    {
        //public Node Node { get; set; }
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}

