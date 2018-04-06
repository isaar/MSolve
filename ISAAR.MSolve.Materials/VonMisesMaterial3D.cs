// --------------------------------------------------------------------------------------------------------------------
// <copyright file="VonMisesMaterial3D.cs" company="National Technical University of Athens">
//   To be decided
// </copyright>
// <summary>
//   Class for 3D Von Mises materials.
// </summary>
// --------------------------------------------------------------------------------------------------------------------
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;


namespace ISAAR.MSolve.FEM.Materials
{
    /// <summary>
    ///   Class for 3D Von Mises materials.
    /// </summary>
    /// <a href = "http://en.wikipedia.org/wiki/Von_Mises_yield_criterion">Wikipedia -Von Mises yield criterion</a>
    public class VonMisesMaterial3D : IIsotropicContinuumMaterial3D
    {
        /// <summary>
        ///   The Poisson ratio value of an incompressible solid.
        /// </summary>
        private const double PoissonRatioForIncompressibleSolid = 0.5;

        /// <summary>
        ///   The total number of strains.
        /// </summary>
        private const int TotalStrains = 6;

        /// <summary>
        ///   The total number of stresses.
        /// </summary>
        private const int TotalStresses = TotalStrains;

        /// <summary>
        ///   An array needed for the formulation of the consistent constitutive matrix.
        /// </summary>
        private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrix = new[,]
            {
                {
                   0, -1, -1, 0, 0, 0
                }, {
                      -1, 0, -1, 0, 0, 0,
                   }, {
                         -1, -1, 0, 0, 0, 0
                      }, {
                            0, 0, 0, 0.5, 0, 0
                         },
                {
                   0, 0, 0, 0, 0.5, 0
                }, {
                      0, 0, 0, 0, 0, 0.5
                   }
            };

        /// <summary>
        ///   The constitutive matrix of the material while still in the elastic region.
        /// </summary>
        private readonly double[,] elasticConstitutiveMatrix;

        /// <summary>
        ///   Hardening modulus for linear hardening.
        /// </summary>
        private readonly double hardeningRatio;

        /// <summary>
        ///   The Poisson ratio.
        /// </summary>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
        /// </remarks>
        private double poissonRatio;

        /// <summary>
        ///   The shear modulus.
        /// </summary>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Shear_modulus">Wikipedia - Shear Modulus</a>
        /// </remarks>
        private readonly double shearModulus;

        /// <summary>
        ///   The yields stress.
        /// </summary>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Yield_%28engineering%29">Yield (engineering)</a>
        ///   The yield strength or yield point of a material is defined in engineering and materials science as the stress at which a material begins to deform plastically.
        /// </remarks>
        private readonly double yieldStress;

        /// <summary>
        ///   The young modulus.
        /// </summary>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
        /// </remarks>
        private double youngModulus;

        /// <summary>
        ///   The constitutive matrix of the material.
        /// </summary>
        private double[,] constitutiveMatrix;

        /// <summary>
        ///   The array of incremental strains.
        /// </summary>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
        /// </remarks>
        private double[] incrementalStrains = new double[6];

        /// <summary>
        ///   Indicates whether this <see cref = "IFiniteElementMaterial" /> is modified.
        /// </summary>
        private bool modified;

        /// <summary>
        ///   The current plastic strain.
        /// </summary>
        private double plasticStrain;

        /// <summary>
        ///   The new plastic strain.
        /// </summary>
        private double plasticStrainNew;

        /// <summary>
        ///   The array of stresses.
        /// </summary>
        private double[] stresses = new double[6];

        /// <summary>
        ///   The array of new stresses.
        /// </summary>
        private double[] stressesNew = new double[6];

        /// <summary>
        ///   Initializes a new instance of the <see cref = "VonMisesMaterial3D" /> class.
        /// </summary>
        /// <param name = "youngModulus">
        ///   The young modulus.
        /// </param>
        /// <param name = "poissonRatio">
        ///   The Poisson ratio.
        /// </param>
        /// <param name = "yieldStress">
        ///   The yield stress.
        /// </param>
        /// <param name = "hardeningRatio">
        ///   The hardening ratio.
        /// </param>
        /// <exception cref = "ArgumentException"> When Poisson ratio is equal to 0.5.</exception>
        public VonMisesMaterial3D(double youngModulus, double poissonRatio, double yieldStress, double hardeningRatio)
        {
            this.youngModulus = youngModulus;

            if (poissonRatio == PoissonRatioForIncompressibleSolid)
            {
                throw new ArgumentException(
                    "Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
            }

            this.poissonRatio = poissonRatio;
            this.yieldStress = yieldStress;
            this.hardeningRatio = hardeningRatio;

            this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
            double lamda = (youngModulus * poissonRatio) / ((1 + poissonRatio) * (1 - (2 * poissonRatio)));
            double mi = youngModulus / (2 * (1 + poissonRatio));
            double value1 = (2 * mi) + lamda;

            this.elasticConstitutiveMatrix = new double[6, 6];
            this.elasticConstitutiveMatrix[0, 0] = value1;
            this.elasticConstitutiveMatrix[0, 1] = lamda;
            this.elasticConstitutiveMatrix[0, 2] = lamda;
            this.elasticConstitutiveMatrix[1, 0] = lamda;
            this.elasticConstitutiveMatrix[1, 1] = value1;
            this.elasticConstitutiveMatrix[1, 2] = lamda;
            this.elasticConstitutiveMatrix[2, 0] = lamda;
            this.elasticConstitutiveMatrix[2, 1] = lamda;
            this.elasticConstitutiveMatrix[2, 2] = value1;
            this.elasticConstitutiveMatrix[3, 3] = mi;
            this.elasticConstitutiveMatrix[4, 4] = mi;
            this.elasticConstitutiveMatrix[5, 5] = mi;
        }

        public double[] Coordinates { get; set; }

        /// <summary>
        ///   Gets the constitutive matrix.
        /// </summary>
        /// <value>
        ///   The constitutive matrix.
        /// </value>
        public ElasticityTensorContinuum3D ConstitutiveMatrix
        {
            get
            {
                if (this.constitutiveMatrix == null)
                {
                    this.UpdateMaterial(new StressStrainVectorContinuum3D(new double[6]));
                }

                return new ElasticityTensorContinuum3D(this.constitutiveMatrix);
            }
        }

        /// <summary>
        ///   Gets the ID of the material.
        /// </summary>
        /// <value>
        ///   The id.
        /// </value>
        public int ID
        {
            get
            {
                return 1;
            }
        }

        /// <summary>
        ///   Gets the incremental strains of the finite element's material.
        /// </summary>
        /// <value>
        ///   The incremental strains.
        /// </value>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
        /// </remarks>
        public double[] IncrementalStrains
        {
            get
            {
                return this.incrementalStrains;
            }
        }

        /// <summary>
        ///   Gets a value indicating whether this <see cref = "IFiniteElementMaterial" /> is modified.
        /// </summary>
        /// <value>
        ///   <c>true</c> if modified; otherwise, <c>false</c>.
        /// </value>
        public bool Modified
        {
            get
            {
                return this.modified;
            }
        }

        /// <summary>
        ///   Gets the plastic strain.
        /// </summary>
        /// <value>
        ///   The plastic strain.
        /// </value>
        public double PlasticStrain
        {
            get
            {
                return this.plasticStrain;
            }
        }

        /// <summary>
        ///   Gets the Poisson ratio.
        /// </summary>
        /// <value>
        ///   The Poisson ratio.
        /// </value>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
        /// </remarks>
        public double PoissonRatio
        {
            get
            {
                return this.poissonRatio;
            }
            set
            {
                this.poissonRatio = value;
            }
        }

        /// <summary>
        ///   Gets the stresses of the finite element's material.
        /// </summary>
        /// <value>
        ///   The stresses.
        /// </value>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Stress_%28mechanics%29">Stress (mechanics)</a>
        /// </remarks>
        public StressStrainVectorContinuum3D Stresses
        {
            get
            {
                return new StressStrainVectorContinuum3D(this.stressesNew);
            }
        }

        /// <summary>
        ///   Gets the Young's Modulus.
        /// </summary>
        /// <value>
        ///   The young modulus.
        /// </value>
        /// <remarks>
        ///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
        /// </remarks>
        public double YoungModulus
        {
            get
            {
                return this.youngModulus;
            }
            set
            {
                this.youngModulus = value;
            }
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        public object Clone()
        {
            var strainsCopy = new double[incrementalStrains.Length];
            Array.Copy(incrementalStrains, strainsCopy, incrementalStrains.Length);
            var stressesCopy = new double[stresses.Length];
            Array.Copy(stresses, stressesCopy, stresses.Length);

            VonMisesMaterial3D m = new VonMisesMaterial3D(
                this.youngModulus, this.poissonRatio, this.yieldStress, this.hardeningRatio)
            {
                modified = this.Modified,
                plasticStrain = this.plasticStrain,
                incrementalStrains = strainsCopy,
                stresses = stressesCopy
            };
            return m;
        }

        /// <summary>
        ///   Resets the indicator of whether the material is modified.
        /// </summary>
        public void ResetModified()
        {
            this.modified = false;
        }

        /// <summary>
        ///   Clears the stresses of the element's material.
        /// </summary>
        public void ClearStresses()
        {
            Array.Clear(this.stresses, 0, 6);
            Array.Clear(this.stressesNew, 0, 6);
        }

        public void ClearState()
        {
            modified = false;
            Array.Clear(constitutiveMatrix, 0, constitutiveMatrix.Length);
            Array.Clear(incrementalStrains, 0, incrementalStrains.Length);
            Array.Clear(stresses, 0, stresses.Length);
            Array.Clear(stressesNew, 0, stressesNew.Length);
            plasticStrain = 0;
            plasticStrainNew = 0;
        }

        /// <summary>
        ///   Saves the state of the element's material.
        /// </summary>
        public void SaveState()
        {
            this.plasticStrain = this.plasticStrainNew;
            Array.Copy(this.stressesNew, this.stresses, 6);
            //this.stresses = this.stressesNew.DeepClone();
        }

        /// <summary>
        ///   Updates the element's material with the provided incremental strains.
        /// </summary>
        /// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
        public void UpdateMaterial(StressStrainVectorContinuum3D strainsIncrement)
        {
            Array.Copy(strainsIncrement.Data, this.incrementalStrains, 6);
            //this.incrementalStrains = strainsIncrement.DeepClone();
            this.CalculateNextStressStrainPoint();
        }

        /// <summary>
        ///   Builds the consistent tangential constitutive matrix.
        /// </summary>
        /// <param name = "value1"> This is a constant already calculated in the calling method. </param>
        /// <remarks>
        ///   We need an additional constant here equal to: 2*G*G*deltaPlasticStrain*sqrt(3/J2elastic)
        ///   Since value1 is already calculated, the additional constant can be calculated through it:
        ///   value3 = 2 * G * value1;
        /// </remarks>
        private void BuildConsistentTangentialConstitutiveMatrix(double value1)
        {
            this.constitutiveMatrix = new double[TotalStresses, TotalStrains];
            double invariantJ2New = this.GetDeviatorSecondStressInvariant(stressesNew);

            double value2 = (3 * this.shearModulus * this.shearModulus) /
                            ((this.hardeningRatio + (3 * this.shearModulus)) * invariantJ2New);

            double value3 = 2 * this.shearModulus * value1;
            double value4 = 0.5 / invariantJ2New;

            var stressDeviator = this.GetStressDeviator(stressesNew);
            for (int k1 = 0; k1 < TotalStresses; k1++)
            {
                for (int k2 = 0; k2 < TotalStresses; k2++)
                {
                    this.constitutiveMatrix[k2, k1] = this.elasticConstitutiveMatrix[k2, k1] -
                                                      (value2 * stressDeviator[k2] * stressDeviator[k1]) -
                                                      (value3 *
                                                       (SupportiveMatrixForConsistentConstitutiveMatrix[k2, k1] -
                                                        (value4 * stressDeviator[k2] * stressDeviator[k1])));
                }
            }
        }

        /// <summary>
        ///   Builds the tangential constitutive matrix.
        /// </summary>
        private void BuildTangentialConstitutiveMatrix()
        {
            this.constitutiveMatrix = new double[TotalStresses, TotalStrains];
            double invariantJ2New = this.GetDeviatorSecondStressInvariant(stressesNew);

            double value2 = (3 * this.shearModulus * this.shearModulus) /
                            ((this.hardeningRatio + (3 * this.shearModulus)) * invariantJ2New);

            var stressDeviator = this.GetStressDeviator(stressesNew);
            for (int k1 = 0; k1 < TotalStresses; k1++)
            {
                for (int k2 = 0; k2 < TotalStresses; k2++)
                {
                    this.constitutiveMatrix[k2, k1] = this.elasticConstitutiveMatrix[k2, k1] -
                                                      (value2 * stressDeviator[k2] * stressDeviator[k1]);
                }
            }
        }

        /// <summary>
        ///   Calculates the next stress-strain point.
        /// </summary>
        /// <exception cref = "InvalidOperationException"> When the new plastic strain is less than the previous one.</exception>
        private void CalculateNextStressStrainPoint()
        {
            var stressesElastic = new double[6];
            for (int i = 0; i < 6; i++)
            {
                stressesElastic[i] = this.stresses[i];
                for (int j = 0; j < 6; j++)
                    stressesElastic[i] += this.elasticConstitutiveMatrix[i, j] * this.incrementalStrains[j];
            }

            double invariantJ2Elastic = this.GetDeviatorSecondStressInvariant(stressesElastic);
            double vonMisesStress = Math.Sqrt(3 * invariantJ2Elastic);
            double vonMisesStressMinusYieldStress = vonMisesStress -
                                                    (this.yieldStress + (this.hardeningRatio * this.plasticStrain));

            bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;

            if (materialIsInElasticRegion)
            {
                this.stressesNew = stressesElastic;
                this.constitutiveMatrix = this.elasticConstitutiveMatrix;
                this.plasticStrainNew = this.plasticStrain;
            }
            else
            {
                double deltaPlasticStrain = vonMisesStressMinusYieldStress /
                                            ((3 * this.shearModulus) + this.hardeningRatio);
                this.plasticStrainNew = this.plasticStrain + deltaPlasticStrain;

                // 2.0 and 1/2 cancel out
                double value1 = this.shearModulus * deltaPlasticStrain * Math.Sqrt(3 / invariantJ2Elastic);
                var stressDeviatorElastic = this.GetStressDeviator(stressesElastic);
                for (int i = 0; i < 6; i++)
                    this.stressesNew[i] = stressesElastic[i] - (value1 * stressDeviatorElastic[i]);

                //this.BuildConsistentTangentialConstitutiveMatrix(value1);
                this.BuildTangentialConstitutiveMatrix();
            }

            if (Math.Abs(this.plasticStrainNew) < Math.Abs(this.plasticStrain))
            {
                throw new InvalidOperationException("Plastic strain cannot decrease.");
            }

            this.modified = this.plasticStrainNew != this.plasticStrain;
        }

        /// <summary>
        ///   Calculates and returns the first stress invariant (I1).
        /// </summary>
        /// <returns> The first stress invariant (I1).</returns>
        public double GetFirstStressInvariant(double[] stresses)
        {
            return stresses[0] + stresses[1] + stresses[2];
        }

        /// <summary>
        ///   Calculates and returns the mean hydrostatic stress.
        /// </summary>
        /// <returns> The mean hydrostatic stress.</returns>
        public double GetMeanStress(double[] stresses)
        {
            double i1 = this.GetFirstStressInvariant(stresses);
            double mean = i1 / 3;
            return mean;
        }

        /// <summary>
        ///   Calculates and returns the second stress invariant (I2).
        /// </summary>
        /// <returns> The second stress invariant (I2).</returns>
        public double GetSecondStressInvariant(double[] stresses)
        {
            return (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2]) - Math.Pow(stresses[5], 2) -
                   Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);
        }

        /// <summary>
        ///   Calculates and returns the stress deviator tensor in vector form.
        /// </summary>
        /// <returns> The stress deviator tensor in vector form.</returns>
        public double[] GetStressDeviator(double[] stresses)
        {
            var hydrostaticStress = this.GetMeanStress(stresses);
            var stressDeviator = new double[]
            {
                stresses[0] - hydrostaticStress,
                stresses[1] - hydrostaticStress,
                stresses[2] - hydrostaticStress,
                stresses[3],
                stresses[4],
                stresses[5]
            };

            return stressDeviator;
        }

        /// <summary>
        ///   Calculates and returns the third stress invariant (I3).
        /// </summary>
        /// <returns> The third stress invariant (I3). </returns>
        public double GetThirdStressInvariant(double[] stresses)
        {
            return (stresses[0] * stresses[1] * stresses[2]) + (2 * stresses[5] * stresses[3] * stresses[4]) - (Math.Pow(stresses[5], 2) * stresses[2]) -
                   (Math.Pow(stresses[3], 2) * stresses[0]) - (Math.Pow(stresses[4], 2) * stresses[1]);
        }

        /// <summary>
        ///   Returns the first stress invariant of the stress deviator tensor (J1), which is zero.
        /// </summary>
        /// <returns> The first stress invariant of the stress deviator tensor (J1). </returns>
        public double GetDeviatorFirstStressInvariant(double[] stresses)
        {
            return 0;
        }

        /// <summary>
        ///   Calculates and returns the second stress invariant of the stress deviator tensor (J2).
        /// </summary>
        /// <returns> The second stress invariant of the stress deviator tensor (J2). </returns>
        public double GetDeviatorSecondStressInvariant(double[] stresses)
        {
            double i1 = this.GetFirstStressInvariant(stresses);
            double i2 = this.GetSecondStressInvariant(stresses);

            double j2 = (1 / 3d * Math.Pow(i1, 2)) - i2;
            return j2;
        }

        /// <summary>
        ///   Calculates and returns the third stress invariant of the stress deviator tensor (J3).
        /// </summary>
        /// <returns> The third deviator stress invariant (J3). </returns>
        public double GetDeviatorThirdStressInvariant(double[] stresses)
        {
            double i1 = this.GetFirstStressInvariant(stresses);
            double i2 = this.GetSecondStressInvariant(stresses);
            double i3 = this.GetThirdStressInvariant(stresses);

            double j3 = (2 / 27 * Math.Pow(i1, 3)) - (1 / 3 * i1 * i2) + i3;
            return j3;
        }

    }
}