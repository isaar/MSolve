/*
 * Copyright © 2006 Stefan Troschütz (stefan@troschuetz.de)
 * 
 * This file is part of Troschuetz.Random Class Library.
 * 
 * Troschuetz.Random is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 * FisherSnedecorDistribution.cs, 27.03.2007
 * 
 * 16.08.2006: Initial version
 * 27.03.2007: Overridden the now virtual base class method Reset to explicitely reset the underlying 
 *				 ChiSquareDistributions of the FisherSnedecorDistribution
 * 
 */

using System;
using Troschuetz.Random;

namespace Troschuetz.Random
{
    /// <summary>
	/// Provides generation of Fisher-Snedecor distributed random numbers.
	/// </summary>
	/// <remarks>
    /// The implementation of the <see cref="FisherSnedecorDistribution"/> type bases upon information presented on
    ///   <a href="http://en.wikipedia.org/wiki/F-distribution">Wikipedia - F-distribution</a>.
    /// </remarks>
	public class FisherSnedecorDistribution : Distribution
	{
		#region instance fields
		/// <summary>
		/// Gets or sets the parameter alpha which is used for generation of Fisher-Snedecor distributed random numbers.
		/// </summary>
		/// <remarks>Call <see cref="IsValidAlpha"/> to determine whether a value is valid and therefor assignable.</remarks>
		public int Alpha
		{
			get
			{
				return this.alpha;
			}
			set
			{
                if (this.IsValidAlpha(value))
                {
                    this.alpha = value;
                    this.UpdateHelpers();
                }
        	}
		}

		/// <summary>
		/// Stores the parameter alpha which is used for generation of Fisher-Snedecor distributed random numbers.
		/// </summary>
        private int alpha;
		
		/// <summary>
        /// Gets or sets the parameter beta which is used for generation of Fisher-Snedecor distributed random numbers.
		/// </summary>
		/// <remarks>Call <see cref="IsValidBeta"/> to determine whether a value is valid and therefor assignable.</remarks>
        public int Beta
		{
			get
			{
				return this.beta;
			}
			set
			{
                if (this.IsValidBeta(value))
                {
                    this.beta = value;
                    this.UpdateHelpers();
                }
        	}
		}

		/// <summary>
        /// Stores the parameter beta which is used for generation of Fisher-Snedecor distributed random numbers.
		/// </summary>
        private int beta;

        /// <summary>
        /// Stores a <see cref="ChiSquareDistribution"/> object used for generation of Fisher-Snedecor distributed random 
        ///   numbers and configured with parameter <see cref="alpha"/>.
        /// </summary>
        private ChiSquareDistribution chiSquareDistributionAlpha;

        /// <summary>
        /// Stores a <see cref="ChiSquareDistribution"/> object used for generation of Fisher-Snedecor distributed random 
        ///   numbers and configured with parameter <see cref="beta"/>.
        /// </summary>
        private ChiSquareDistribution chiSquareDistributionBeta;

        /// <summary>
        /// Stores an intermediate result for generation of Fisher-Snedecor distributed random numbers.
        /// </summary>
        /// <remarks>
        /// Speeds up random number generation cause this value only depends on distribution parameters 
        ///   and therefor doesn't need to be recalculated in successive executions of <see cref="NextDouble"/>.
        /// </remarks>
        private double helper1;
        #endregion

		#region construction
		/// <summary>
        /// Initializes a new instance of the <see cref="FisherSnedecorDistribution"/> class, using a 
        ///   <see cref="StandardGenerator"/> as underlying random number generator.
		/// </summary>
        public FisherSnedecorDistribution()
            : this(new StandardGenerator())
		{
		}

		/// <summary>
        /// Initializes a new instance of the <see cref="FisherSnedecorDistribution"/> class, using the specified 
        ///   <see cref="Generator"/> as underlying random number generator.
        /// </summary>
        /// <param name="generator">A <see cref="Generator"/> object.</param>
        /// <exception cref="ArgumentNullException">
        /// <paramref name="generator"/> is NULL (<see langword="Nothing"/> in Visual Basic).
        /// </exception>
        public FisherSnedecorDistribution(Generator generator)
            : base(generator)
        {
            this.alpha = 1;
            this.beta = 1;
            this.chiSquareDistributionAlpha = new ChiSquareDistribution(generator);
            this.chiSquareDistributionBeta = new ChiSquareDistribution(generator);
            this.UpdateHelpers();
        }
		#endregion
	
		#region instance methods
		/// <summary>
		/// Determines whether the specified value is valid for parameter <see cref="Alpha"/>.
		/// </summary>
		/// <param name="value">The value to check.</param>
		/// <returns>
		/// <see langword="true"/> if value is greater than 0; otherwise, <see langword="false"/>.
		/// </returns>
        public bool IsValidAlpha(int value)
		{
			return value > 0;
		}
		
		/// <summary>
		/// Determines whether the specified value is valid for parameter <see cref="Beta"/>.
		/// </summary>
		/// <param name="value">The value to check.</param>
		/// <returns>
		/// <see langword="true"/> if value is greater than 0; otherwise, <see langword="false"/>.
		/// </returns>
        public bool IsValidBeta(int value)
		{
			return value > 0;
		}
        
        /// <summary>
        /// Updates the helper variables that store intermediate results for generation of Fisher-Snedecor distributed random 
        ///   numbers.
        /// </summary>
        private void UpdateHelpers()
        {
            this.chiSquareDistributionAlpha.Alpha = this.alpha;
            this.chiSquareDistributionBeta.Alpha = this.beta;
            this.helper1 = (double)this.beta / (double)this.alpha;
        }
        #endregion

		#region overridden Distribution members
		/// <summary>
		/// Resets the Fisher-Snedecor distribution, so that it produces the same random number sequence again.
		/// </summary>
		/// <returns>
		/// <see langword="true"/>, if the Fisher-Snedecor distribution was reset; otherwise, <see langword="false"/>.
		/// </returns>
		public override bool Reset()
		{
			bool result = base.Reset();
			if (result)
			{
				result = this.chiSquareDistributionAlpha.Reset();
				if (result)
				{
					result = this.chiSquareDistributionBeta.Reset();
				}
			}

			return result;
		}

		/// <summary>
		/// Gets the minimum possible value of Fisher-Snedecor distributed random numbers.
		/// </summary>
		public override double Minimum
		{
			get
			{
				return 0.0;
			}
		}

		/// <summary>
		/// Gets the maximum possible value of Fisher-Snedecor distributed random numbers.
		/// </summary>
        public override double Maximum
		{
			get
			{
				return double.MaxValue;
			}
		}

		/// <summary>
		/// Gets the mean value of Fisher-Snedecor distributed random numbers.
		/// </summary>
        public override double Mean
		{
			get
			{
                if (this.beta > 2)
                {
                    return (double)this.beta / ((double)this.beta - 2.0);
                }
                else
                {
                    return double.NaN;
                }
			}
		}
		
		/// <summary>
		/// Gets the median of Fisher-Snedecor distributed random numbers.
		/// </summary>
        public override double Median
		{
			get
			{
                return double.NaN;
			}
		}
		
		/// <summary>
		/// Gets the variance of Fisher-Snedecor distributed random numbers.
		/// </summary>
        public override double Variance
		{
			get
			{
                if (this.beta > 4)
                {
                    double alpha = this.alpha;
                    double beta = this.beta;

                    return 2 * Math.Pow(this.beta, 2.0) *(alpha + beta - 2.0) / alpha / Math.Pow(this.beta - 2.0, 2.0) / 
                        (this.beta - 4.0);
                }
                else
                {
                    return double.NaN;
                }
			}
		}
		
		/// <summary>
		/// Gets the mode of Fisher-Snedecor distributed random numbers.
		/// </summary>
        public override double[] Mode
		{
			get
			{
                if (this.alpha > 2)
                {
                    double alpha = this.alpha;
                    double beta = this.beta;

                    return new double[] { (alpha - 2.0) / alpha * beta / (beta + 2.0) };
                }
                else
                {
                    return new double[] { };
                }
			}
		}
		
		/// <summary>
		/// Returns a Fisher-Snedecor distributed floating point random number.
		/// </summary>
		/// <returns>A Fisher-Snedecor distributed double-precision floating point number.</returns>
        public override double NextDouble()
		{
            return this.chiSquareDistributionAlpha.NextDouble() / this.chiSquareDistributionBeta.NextDouble() * this.helper1;
		}
		#endregion
	}
}