using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ThermalMaterial
    {
        public ThermalMaterial(double density, double specialHeatCoeff, double thermalConductivity)
        {
            this.Density = density;
            this.SpecialHeatCoeff = specialHeatCoeff;
            this.ThermalConductivity = thermalConductivity;
        }

        public double Density { get; }
        public double SpecialHeatCoeff { get; }
        public double ThermalConductivity { get; }

        public ThermalMaterial Clone() => new ThermalMaterial(Density, SpecialHeatCoeff, ThermalConductivity);
    }
}
