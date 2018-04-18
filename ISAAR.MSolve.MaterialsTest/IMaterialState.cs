using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.MaterialsTest
{
    public interface IMaterialState
    {
        bool Modified { get; }
        void ResetModified();
        double[] Stresses { get; }
        double[] ConstitutiveMatrix { get; }
        void UpdateMaterial(double[] strains);
        void ClearState();
        void SaveState();
        void ClearStresses();
    }
}
