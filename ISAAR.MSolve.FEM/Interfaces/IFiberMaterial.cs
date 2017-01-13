using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IFiberMaterial : IFiniteElementMaterial
    {
        double Stress { get; }
        double Strain { get; }
        void UpdateMaterial(double dStrain);
        void SaveState();
        void ClearStresses();
        IFiberMaterial Clone(IFiberFiniteElementMaterial parent);
        //double YoungModulusElastic { get; set; }
    }
}
