using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Entities
{
    public class MassAccelerationSinusoidalLoad : IMassAccelerationHistoryLoad
    {
        private readonly double period, magnitude, timeStep, timeLength;
        private readonly int timeSteps;
        public virtual DOFType DOF { get; set; }

        public virtual double this[int currentTimeStep]
        {
            get { return currentTimeStep < timeSteps ? magnitude * Math.Sin(2 * Math.PI / period * timeStep * currentTimeStep) : 0; }
        }

        public MassAccelerationSinusoidalLoad(double period, double timeStep, double timeLength, double magnitude)
        {
            this.period = period;
            this.timeStep = timeStep;
            this.timeLength = timeLength;
            this.magnitude = magnitude;
            this.timeSteps = (int)(timeLength / timeStep);
        }

    }

}
