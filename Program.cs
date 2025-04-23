using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Mime;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using OpenTDv241;
using static OpenTDv241.Results.Plot.Axis;


namespace VectorGenerator
{
    public class Program
    {
        double h = 500;
        double beta_deg = 10 * -1;
        double beta_rad;
        List<double> times = new List<double>();
        List<Vector3d> positions = new List<Vector3d>();
        List<Vector3d> velocities = new List<Vector3d>();
        List<int> attFlag = new List<int>();
        List<Vector3d> GCS_nadir_points = new List<Vector3d>();
        Vector3d GCS_sun_point = new Vector3d(0, 1, 0);
        Vector3d GCS_heavens_point = new Vector3d(0, 0, 1);
        Vector3d x0 = new Vector3d(1, 0, 0);
        Vector3d y0 = new Vector3d(0, 1, 0);
        Vector3d z0 = new Vector3d(0, 0, 1);
        List<Vector3d> GCS_basis_matrix = new List<Vector3d>();
        List<Vector3d> planet_vectors = new List<Vector3d>();
        List<Vector3d> sun_vectors = new List<Vector3d>();
        List<List<Vector3d>> SCCS_basis_matrices = new List<List<Vector3d>>();
        List<double> radii = new List<double>();

        Vector3d V_no;

        double R = 6378.0;
        double mu = 398600;
        double a;
        double T_orbit_sec;
        double beta_star;

        double lastTime = 0;
        double eclipse_fraction = 0;
        double t_shadow_entry = 0;
        double t_shadow_exit = 0;
        double eclipse_duration = 0;

        public void Run()
        {
            a = R + h;
            T_orbit_sec = Math.Floor(2 * Math.PI * Math.Pow(mu, -0.5) * Math.Pow(a, 1.5));
            beta_star = Math.Asin(R / (R + h));
            beta_rad = beta_deg * Math.PI / 180.0;
            //times.Add(0.0);
            //attFlag.Add(1); 

            // orbit normal vector
            V_no = new Vector3d( 0, -1 * Math.Sin(beta_rad), Math.Cos(beta_rad) );
            
            // global coordinate system basis matrix
            GCS_basis_matrix.Add(x0);
            GCS_basis_matrix.Add(y0);
            GCS_basis_matrix.Add(z0);

            
            // updating eclipse times, use variables as needed when defining attitude durations
            if (beta_rad < beta_star)
            {
                eclipse_fraction = Math.Acos(Math.Sqrt(Math.Pow(h, 2) + 2 * R * h) / ((R + h) * Math.Cos(beta_rad))) / Math.PI;
                eclipse_duration = eclipse_fraction * T_orbit_sec;
                t_shadow_entry = T_orbit_sec / 2 - eclipse_duration / 2;
                t_shadow_exit = t_shadow_entry + eclipse_duration;
            }

            
            // add attitudes here !!
            addSunpoint(T_orbit_sec / 3, 5);
            addAuxiliarypoint(T_orbit_sec / 3, 5, 20);
            addSunpoint(T_orbit_sec / 3, 5);
            calc_positions();
            calc_pointing_vectors();
            calc_velocity_vectors();
            calc_SC_basis_matrices(); 
            calc_output_arrays();
            debugPrints();

            // replace path variable with dwg path from your file system
            var path = @"C:\Users\AndrewDevine\Documents\Models\GEO_Aries\safe - Copy.dwg";
            // name the orbit what you'd like
            var orbitName = "Trajectory_custom";
            var td = new ThermalDesktop(path);
            td.Connect();
            var testOrbit = td.CreateOrbit(orbitName);
            testOrbit.OrbitType = OpenTDv241.RadCAD.Orbit.OrbitTypes.TRAJECTORY;
            testOrbit.HrTimeArray = times;
            testOrbit.HrPlanetVecArray = planet_vectors;
            testOrbit.HrSunVecArray = sun_vectors;
            testOrbit.HrOrbitRadiusArray = radii;
            testOrbit.Update();
            Console.WriteLine("program finished, press enter to end");
            Console.ReadLine();
        }

        public static void Main(string[] args)
        {
            new Program().Run();
        }

        /// <summary>
        /// creates output arrays to be inserted into TD by multiplying the SC basis matrices by global pointing vectors
        /// </summary>
        public void calc_output_arrays()
        {

            for (var i = 0; i < positions.Count; i++)
            {
                sun_vectors.Add(MultiplyMatrixByVector(SCCS_basis_matrices[i], GCS_sun_point));
                planet_vectors.Add(MultiplyMatrixByVector(SCCS_basis_matrices[i], GCS_nadir_points[i]));
            }
        }

        /// <summary>
        /// creates SC basis matrices by aligning coordinate axes to global vectors
        /// </summary>
        public void calc_SC_basis_matrices()
        {
            for (var i = 0; i < positions.Count; i++)
            {
                if (attFlag[i] == 1)
                {
                    SCCS_basis_matrices[i][2] = -1 * GCS_sun_point;
                    SCCS_basis_matrices[i][1] = -1 * GCS_heavens_point;
                    SCCS_basis_matrices[i][0] = SCCS_basis_matrices[i][2].CrossProduct(SCCS_basis_matrices[i][1]);
                } 
                else if (attFlag[i] == 2)
                {
                    var z = GCS_nadir_points[i];
                    var x = velocities[i];
                    var y = -1 * x.CrossProduct(z);

                    x.Normalize();
                    z.Normalize();
                    y.Normalize();
                    SCCS_basis_matrices[i][0] = x;
                    SCCS_basis_matrices[i][1] = y;
                    SCCS_basis_matrices[i][2] = z;
                }
                else
                {
                    var angleRad = attFlag[i] * Math.PI / 180.0;
                    SCCS_basis_matrices[i][0] = velocities[i];
                    var z = GCS_nadir_points[i];
                    var x = velocities[i];
                    var z_rotated = z * Math.Cos(angleRad) + x.CrossProduct(z) * Math.Sin(angleRad) + x * x.DotProduct(z) * (1 - Math.Cos(angleRad));
                    SCCS_basis_matrices[i][2] = z_rotated;
                    SCCS_basis_matrices[i][1] = 1 * velocities[i].CrossProduct(z_rotated);
                }
            }
        }

        /// <summary>
        /// calculates velocity vectors by cross product of position and orbit normal vectors
        /// </summary>
        public void calc_velocity_vectors()
        {
            foreach (Vector3d position in positions)
            {
                var velocity = V_no.CrossProduct(position);
                velocity.Normalize(); 
                velocities.Add(velocity);
            }
        }
        /// <summary>
        /// calculates 3D positions, fills radii array and initializes SC basis matrices to zero
        /// </summary>
        public void calc_positions()
        {
            foreach (double time in times)
            {
                positions.Add(new Vector3d( getX(time), getY(time), getZ(time)));
                radii.Add(a / R);
                var basis_matrix = new List<Vector3d>();

                for (var j = 0; j < 3; j++)
                {

                    Vector3d vector = new Vector3d(0, 0, 0);
                    basis_matrix.Add(vector);

                }

                SCCS_basis_matrices.Add(basis_matrix);
            }
           
        }
        /// <summary>
        /// creates nadir and solar pointing vectors in global coordinates
        /// </summary>
        public void calc_pointing_vectors() 
        {
            foreach (Vector3d position in positions) 
            {
                var nadirpoint = new Vector3d(0,0,0) - position;
                nadirpoint.Normalize();
                GCS_nadir_points.Add(nadirpoint);
            }
        }
        /// <summary>
        /// parametric x position
        /// </summary>
        /// <param name="t">time, in seconds</param>
        /// <returns> x position at time t</returns>
        public double getX(double t) 
        {

            return (a * Math.Cos(t / T_orbit_sec * 2 * Math.PI + Math.PI / 2));
       
        }
        /// <summary>
        /// parametric y position
        /// </summary>
        /// <param name="t">time, in seconds</param>
        /// <returns> y position at time t</returns>
        public double getY(double t)
        {

            return (a * Math.Sin(t / T_orbit_sec * 2 * Math.PI + Math.PI / 2) * Math.Cos(beta_rad));
         
        }
        /// <summary>
        /// parametric z position
        /// </summary>
        /// <param name="t">time, in seconds</param>
        /// <returns> z position at time t</returns>
        public double getZ(double t)
        {

            return (a * Math.Sin(t / T_orbit_sec * 2 * Math.PI + Math.PI / 2) * Math.Sin(beta_rad));

        }

        /// <summary>
        /// this method adds a sun-pointing attitude to the orbit. 
        /// the solar arrays (Z- of SC) will be aligned with the sun.
        /// the SC Y- will be perpendicular to the ecliptic plane.
        /// </summary>
        /// <param name="duration">duration of the attitude, in seconds</param>
        /// <param name="resolution">quantity of positions for the attitude, use an integer</param>
        public void addSunpoint(double duration, int resolution)
        {
            resolution++;
            double currentTime = lastTime + duration;
            double[] spTimes = Generate.LinearSpaced(resolution, lastTime, currentTime);
            lastTime = currentTime;
            times.AddRange(spTimes.Skip(1));
            attFlag.AddRange(Enumerable.Repeat(1, resolution-1));
        }

        /// <summary>
        ///  this method adds a nadir-pointing attitude to the orbit.
        ///  the X+ of the SC will be constrained to velocity.
        ///  the Z+ of the SC will be pointing nadir.
        /// </summary>
        /// <param name="duration">duration of the attitude, in seconds </param>
        /// <param name="resolution">quantity of positions for the attitude, use an integer</param>
        public void addNadirpoint(double duration, int resolution)
        {
            resolution++;
            double currentTime = lastTime + duration;
            double[] npTimes = Generate.LinearSpaced(resolution, lastTime, currentTime);
            lastTime = currentTime;
            times.AddRange(npTimes.Skip(1));
            attFlag.AddRange(Enumerable.Repeat(2, resolution-1));
        }

        /// <summary>
        /// this method adds an off-pointing attitude to the orbit.
        /// this off-pointing angle is calculated by using the cross product of the nadir vector and the velocity vector.
        /// the X+ of the spacecraft will be constrained to velocity.
        /// the Z+ of the spacecraft will be constrained to nadir + the input angle.
        /// </summary>
        /// <param name="duration">duration of the attitude, in second</param>
        /// <param name="resolution">quantity of positions for the attitude, use an integer</param>
        /// <param name="angleDeg">off-angle of the attitude, in degrees</param>
        public void addAuxiliarypoint(double duration, int resolution, int angleDeg)
        {
            resolution++;
            double currentTime = lastTime + duration;
            double[] apTimes = Generate.LinearSpaced(resolution, lastTime, currentTime);
            lastTime = currentTime;
            times.AddRange(apTimes.Skip(1));
            attFlag.AddRange(Enumerable.Repeat(angleDeg, resolution-1));
        }


        // internal math, transposes and multiplies
        public static Vector3d MultiplyMatrixByVector(List<Vector3d> matrix, Vector3d vector)
        {
            // Treat the matrix as row vectors to convert from global to local frame (transpose)
            double x = vector.DotProduct(matrix[0]); // row 0 = x-axis of SCCS in GCS
            double y = vector.DotProduct(matrix[1]); // row 1 = y-axis of SCCS in GCS
            double z = vector.DotProduct(matrix[2]); // row 2 = z-axis of SCCS in GCS

            return new Vector3d(x, y, z);
        }

        // do not use, wrong dimensions
        //public static Vector3d MultiplyMatrixByVectorSolar(List<Vector3d> matrix, Vector3d vector)
        //{
        //    if (matrix.Count != 3)
        //        throw new ArgumentException("Matrix must have exactly 3 Vector3d columns.");

        //    double x =
        //        matrix[0].X * vector.X +
        //        matrix[1].X * vector.Y +
        //        matrix[2].X * vector.Z;

        //    double y =
        //        matrix[0].Y * vector.X +
        //        matrix[1].Y * vector.Y +
        //        matrix[2].Y * vector.Z;

        //    double z =
        //        matrix[0].Z * vector.X +
        //        matrix[1].Z * vector.Y +
        //        matrix[2].Z * vector.Z;

        //    return new Vector3d(x, y, z);
        //}

        public void debugPrints()
        {
            Console.WriteLine("Orbit period: " + T_orbit_sec);
            Console.WriteLine();
            Console.WriteLine("last time in times: " + times[times.Count - 1]);
            Console.WriteLine("");
            Console.WriteLine("times: ");

            for (var i= 0; i < times.Count; i++)
            {
                Console.WriteLine("time: " + times[i] + ", flag: " + attFlag[i]); 
            }
            Console.WriteLine("  ");

            
            Console.WriteLine("positions");
            foreach (var position in positions)
            {
                Console.WriteLine(position);
            }
            Console.WriteLine("  ");

            
            Console.WriteLine("nadir points GCS");
            foreach (var vec in GCS_nadir_points)
            {
                Console.WriteLine(vec);
            }
            Console.WriteLine("  ");

            Console.WriteLine("velocities GCS");
            foreach (var vec in velocities)
            {
                Console.WriteLine(vec);
            }
            Console.WriteLine("  ");



            Console.WriteLine(SCCS_basis_matrices);

            Console.WriteLine("SCCS_basis_matrices:  ");
            foreach (var matrix in SCCS_basis_matrices)
            {
                foreach (var vector in matrix)
                {
                    Console.WriteLine(vector);
                }
                Console.WriteLine(' ');

            }

            Console.WriteLine();

            Console.WriteLine("planet_vectors:  ");
            foreach (var vector in planet_vectors)
            {
                Console.WriteLine(vector);

            }

            Console.WriteLine();

            Console.WriteLine("sun_vectors:  ");
            foreach (var vector in sun_vectors)
            {
                Console.WriteLine(vector);

            }
            Console.WriteLine("Length of times: " + times.Count());
            Console.WriteLine("Length of planet vectors: " + planet_vectors.Count());
            Console.WriteLine("Length of solar vectors: " + sun_vectors.Count());
            Console.WriteLine("Length of radiis: " + radii.Count());

        }


    }
}
