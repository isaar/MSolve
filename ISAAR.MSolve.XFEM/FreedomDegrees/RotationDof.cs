using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.FreedomDegrees
{
    public class RotationDof: IDof
    {
        public static readonly RotationDof X = new RotationDof("Rotation around X axis");
        public static readonly RotationDof Y = new RotationDof("Rotation around Y axis");
        public static readonly RotationDof Z = new RotationDof("Rotation around Z axis");

        private readonly string name;

        private RotationDof(string name) // No more displacement dofs can be created
        {
            this.name = name;
        }

        public override string ToString()
        {
            return name;
        }
    }
}
