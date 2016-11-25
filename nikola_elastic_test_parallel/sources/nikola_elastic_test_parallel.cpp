/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004
 */

// @sect3{Include files}

// First the usual assortment of header files we have already used in previous
// example programs:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// And here come the things that we need particularly for this example program
// and that weren't in step-8. First, we replace the standard output
// <code>std::cout</code> by a new stream <code>pcout</code> which is used in
// parallel computations for generating output only on one of the MPI
// processes.
#include <deal.II/base/conditional_ostream.h>
// We are going to query the number of processes and the number of the present
// process by calling the respective functions in the Utilities::MPI
// namespace.
#include <deal.II/base/utilities.h>
// Then, we are going to replace all linear algebra components that involve
// the (global) linear system by classes that wrap interfaces similar to our
// own linear algebra classes around what PETSc offers (PETSc is a library
// written in C, and deal.II comes with wrapper classes that provide the PETSc
// functionality with an interface that is similar to the interface we already
// had for our own linear algebra classes). In particular, we need vectors and
// matrices that are distributed across several processes in MPI programs (and
// simply map to sequential, local vectors and matrices if there is only a
// single process, i.e. if you are running on only one machine, and without
// MPI support):
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
// Then we also need interfaces for solvers and preconditioners that PETSc
// provides:
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
// And in addition, we need some algorithms for partitioning our meshes so
// that they can be efficiently distributed across an MPI network. The
// partitioning algorithm is implemented in the <code>GridTools</code> class,
// and we need an additional include file for a function in
// <code>DoFRenumbering</code> that allows to sort the indices associated with
// degrees of freedom so that they are numbered according to the subdomain
// they are associated with:
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

// And this is simply C++ again:
#include <fstream>
#include <iostream>
#include <sstream>


#include "../../feminist/prmt_sintactic_addition/prmt_sintactic_addition.h"
#include "../../feminist/calculation_core/src/blocks/general/additional_tools/types/types.h"
#include "../../feminist/calculation_core/src/blocks/special/elastic_problem_tools/elastic_problem_tools.h"
#include "../../feminist/calculation_core/src/blocks/special/nikola_problem/source/scalar/source_scalar.h"
#include "../../feminist/calculation_core/src/blocks/special/nikola_problem/source/vector/source_vector.h"
#include "../../feminist/calculation_core/src/blocks/special/feature_source/scalar/source_scalar.h"
#include "../../feminist/calculation_core/src/blocks/special/feature_source/vector/source_vector.h"
#include "../../feminist/calculation_core/src/blocks/special/poly_materials_source/scalar/source_scalar.h"
#include "../../feminist/calculation_core/src/blocks/special/poly_materials_source/vector/source_vector.h"

void set_laminat(dealii::Triangulation< 2 > &triangulation,
        const std::function<dbl(dbl)> delimiter, cdbl lx, cdbl ly, cst np)
{
    enum {x, y, z};

    {
        std::vector< dealii::Point< 2 > > v (8);

        // v[0]  = dealii::Point<2>(0.0, 0.0);
        // v[1]  = dealii::Point<2>(1.0, 0.0);
        // v[2]  = dealii::Point<2>(9.0, 0.0);
        // v[3]  = dealii::Point<2>(10.0, 0.0);
        // v[4]  = dealii::Point<2>(0.0, 1.0);
        // v[5]  = dealii::Point<2>(1.0, 1.0);
        // v[6]  = dealii::Point<2>(9.0, 1.0);
        // v[7]  = dealii::Point<2>(10.0, 1.0);
        v[0]  = dealii::Point<2>(0.0, 0.0);
        v[1]  = dealii::Point<2>(3.3333, 0.0);
        v[2]  = dealii::Point<2>(6.6666, 0.0);
        v[3]  = dealii::Point<2>(10.0, 0.0);
        v[4]  = dealii::Point<2>(0.0, 1.0);
        v[5]  = dealii::Point<2>(3.3333, 1.0);
        v[6]  = dealii::Point<2>(6.6666, 1.0);
        v[7]  = dealii::Point<2>(10.0, 1.0);

        std::vector< dealii::CellData<2>> c; //(3, dealii::CellData<2>());
        {
            dealii::CellData<2> tmp;
            tmp.vertices[0]=0;
            tmp.vertices[1]=1;
            tmp.vertices[2]=4;
            tmp.vertices[3]=5;
            tmp.material_id=0;
            c .push_back (tmp);
            tmp.vertices[0]=1;
            tmp.vertices[1]=2;
            tmp.vertices[2]=5;
            tmp.vertices[3]=6;
            tmp.material_id=0;
            c .push_back (tmp);
            tmp.vertices[0]=2;
            tmp.vertices[1]=3;
            tmp.vertices[2]=6;
            tmp.vertices[3]=7;
            tmp.material_id=0;
            c .push_back (tmp);
        };
        triangulation .create_triangulation (v, c, dealii::SubCellData());
    };
    triangulation.refine_global(9);

    {
        auto cell = triangulation .begin_active();
        auto end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            auto c = cell->center();

            if (c(y) > delimiter(c(x)))
            {
                cell->set_material_id(1);
            }
            else
            {
                cell->set_material_id(0);
            };
        };
    };
};

void get_nikola_stress (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const vec<ATools::FourthOrderTensor> &C,
        arr<dealii::Vector<dbl>, 2> &stress)
{
    enum {x, y, z};

    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        /* Точки 3 и 2 переставленны местами, потому что в диле у них
         * порядок зигзагом, а мне надо по кругу
         */
        arr<prmt::Point<2>, 4> points = {
            prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
            prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
            prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
            prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
        for (st i = 0; i < 4; ++i)
        {
            arr<arr<dbl, 2>, 2> grad;
            for (st component = 0; component < 2; ++component)
            {
                arr<dbl, 4> values = {
                    move(cell->vertex_dof_index (0, component)),
                    move(cell->vertex_dof_index (1, component)),
                    move(cell->vertex_dof_index (3, component)),
                    move(cell->vertex_dof_index (2, component))};

                Scalar4PointsFunc<2> function_on_cell(points, values);
                auto mat_id = cell->material_id();

                grad[0][component] = function_on_cell.dx(cell->vertex(i));
                grad[1][component] = function_on_cell.dy(cell->vertex(i));
            };
            auto indx_1 = cell->vertex_dof_index(i, 0);
            auto indx_2 = cell->vertex_dof_index(i, 1);
            auto mat_id = cell->material_id();

            stress[0][indx_1] += 
                // grad[x][x];//C[mat_id][0][0][0][0];
                C[mat_id][x][x][x][x] * grad[x][x] +
                C[mat_id][x][x][x][y] * grad[x][y] +
                C[mat_id][x][x][y][x] * grad[y][x] +
                C[mat_id][x][x][y][y] * grad[y][y];// + C[mat_id][x][x][z][z];
            stress[x][indx_2] += 
                C[mat_id][x][y][x][x] * grad[x][x] +
                C[mat_id][x][y][x][y] * grad[x][y] +
                C[mat_id][x][y][y][x] * grad[y][x] +
                C[mat_id][x][y][y][y] * grad[y][y];
            stress[y][indx_1] += 
                C[mat_id][y][x][x][x] * grad[x][x] +
                C[mat_id][y][x][x][y] * grad[x][y] +
                C[mat_id][y][x][y][x] * grad[y][x] +
                C[mat_id][y][x][y][y] * grad[y][y];
            stress[y][indx_2] += 
                // grad[y][y];//C[mat_id][0][0][0][0];
                C[mat_id][y][y][x][x] * grad[x][x] +
                C[mat_id][y][y][x][y] * grad[x][y] +
                C[mat_id][y][y][y][x] * grad[y][x] +
                C[mat_id][y][y][y][y] * grad[y][y] + C[mat_id][y][y][z][z];

                ++(N[indx_1]); 
                ++(N[indx_2]); 
            };
    };
    for (st i = 0; i < N.size(); ++i)
    {
        stress[0][i] /= N[i];
        stress[1][i] /= N[i];
    };
};

// The last step is as in all previous programs:
namespace Step17
{
    using namespace dealii;

    template <int dim>
        class ElasticProblem
        {
            public:
                ElasticProblem ();
                ~ElasticProblem ();
                void run ();

            private:
                void setup_system ();
                void assemble_system ();
                unsigned int solve ();
                void refine_grid ();
                void output_results (const unsigned int cycle) const;

                ConditionalOStream pcout;

                Triangulation<dim>   triangulation;
                DoFHandler<dim>      dof_handler;

                FESystem<dim>        fe;

                ConstraintMatrix     hanging_node_constraints;

                typename PETScWrappers::MPI::SparseMatrix system_matrix;

                PETScWrappers::MPI::Vector       solution;
                PETScWrappers::MPI::Vector       system_rhs;

                MPI_Comm mpi_communicator;

                const unsigned int n_mpi_processes;
                const unsigned int this_mpi_process;

                // cdbl E_1   = 1.0;  cdbl E_2   = 10.0;
                // // cdbl E_1   = 1.0;  cdbl E_2   = 1.0;
                // cdbl pua_1 = 0.25; cdbl pua_2 = 0.25;
                cdbl E_1   = 1.0;  cdbl E_2   = 10.0;
                // cdbl E_1   = 1.0;  cdbl E_2   = 1.0;
                cdbl pua_1 = 0.2; cdbl pua_2 = 0.1;
        };

    template <int dim>
        class RightHandSide :  public Function<dim>
    {
        public:
            RightHandSide ();

            virtual void vector_value (const Point<dim> &p,
                    Vector<double>   &values) const;

            virtual void vector_value_list (const std::vector<Point<dim> > &points,
                    std::vector<Vector<double> >   &value_list) const;
    };


    template <int dim>
        RightHandSide<dim>::RightHandSide () :
            Function<dim> (dim)
    {}


    template <int dim>
        inline
        void RightHandSide<dim>::vector_value (const Point<dim> &p,
                Vector<double>   &values) const
        {
            Assert (values.size() == dim,
                    ExcDimensionMismatch (values.size(), dim));
            Assert (dim >= 2, ExcInternalError());

            Point<dim> point_1, point_2;
            point_1(0) = 0.5;
            point_2(0) = -0.5;

            if (((p-point_1).norm_square() < 0.2*0.2) ||
                    ((p-point_2).norm_square() < 0.2*0.2))
                values(0) = 1;
            else
                values(0) = 0;

            if (p.square() < 0.2*0.2)
                values(1) = 1;
            else
                values(1) = 0;
        }

    template <int dim>
        void RightHandSide<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                std::vector<Vector<double> >   &value_list) const
        {
            const unsigned int n_points = points.size();

            Assert (value_list.size() == n_points,
                    ExcDimensionMismatch (value_list.size(), n_points));

            for (unsigned int p=0; p<n_points; ++p)
                RightHandSide<dim>::vector_value (points[p],
                        value_list[p]);
        }



    template <int dim>
        class LambdaValue :  public Function<dim>
    {
        public:
            LambdaValue (cdbl young1, cdbl poisson1, cdbl young2, cdbl poisson2);

            virtual double value (const Point<dim> &p) const;

            virtual void value_list (const std::vector<Point<dim> > &points,
                    std::vector<double> &value_list) const;
            
            dbl lambda1;
            dbl lambda2;
    };


    template <int dim>
        LambdaValue<dim>::LambdaValue (cdbl young1, cdbl poisson1, cdbl young2, cdbl poisson2) :
            Function<dim> (dim)
    {
        lambda1 = (poisson1 * young1) / ((1 + poisson1) * (1 - 2 * poisson1));
        lambda2 = (poisson2 * young2) / ((1 + poisson2) * (1 - 2 * poisson2));
    }


    template <int dim>
        inline
        double LambdaValue<dim>::value (const Point<dim> &p) const
        {
            return p(1) > 0.5 ? lambda2 : lambda1;
        }

    template <int dim>
        void LambdaValue<dim>::value_list (const std::vector<Point<dim> > &points,
                std::vector<double> &value_list) const
        {
            const unsigned int n_points = points.size();

            for (unsigned int p=0; p<n_points; ++p)
                value_list[p] = LambdaValue<dim>::value (points[p]);
        }



    template <int dim>
        class MuValue :  public Function<dim>
    {
        public:
            MuValue (cdbl young1, cdbl poisson1, cdbl young2, cdbl poisson2);

            virtual double value (const Point<dim> &p) const;

            virtual void value_list (const std::vector<Point<dim> > &points,
                    std::vector<double> &value_list) const;
            
            dbl mu1;
            dbl mu2;
    };


    template <int dim>
        MuValue<dim>::MuValue (cdbl young1, cdbl poisson1, cdbl young2, cdbl poisson2) :
            Function<dim> (dim)
    {
        mu1 = young1 / (2 * (1 + poisson1));
        mu2 = young2 / (2 * (1 + poisson2));
    }


    template <int dim>
        inline
        double MuValue<dim>::value (const Point<dim> &p) const
        {
            return p(1) > 0.5 ? mu2 : mu1;
        }

    template <int dim>
        void MuValue<dim>::value_list (const std::vector<Point<dim> > &points,
                std::vector<double> &value_list) const
        {
            const unsigned int n_points = points.size();

            for (unsigned int p=0; p<n_points; ++p)
                value_list[p] = MuValue<dim>::value (points[p]);
        }


    template <int dim>
        ElasticProblem<dim>::ElasticProblem ()
        :
            pcout (std::cout),
            dof_handler (triangulation),
            fe (FE_Q<dim>(1), dim),
            mpi_communicator (MPI_COMM_WORLD),
            n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
            this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator))
    {
        pcout.set_condition(this_mpi_process == 0);
    }


    template <int dim>
        ElasticProblem<dim>::~ElasticProblem ()
        {
            dof_handler.clear ();
        }

    template <int dim>
        void ElasticProblem<dim>::setup_system ()
        {
            GridTools::partition_triangulation (n_mpi_processes, triangulation);

            dof_handler.distribute_dofs (fe);
            DoFRenumbering::subdomain_wise (dof_handler);

            const types::global_dof_index n_local_dofs
                = DoFTools::count_dofs_with_subdomain_association (dof_handler,
                        this_mpi_process);
            system_matrix.reinit (mpi_communicator,
                    dof_handler.n_dofs(),
                    dof_handler.n_dofs(),
                    n_local_dofs,
                    n_local_dofs,
                    dof_handler.max_couplings_between_dofs());

            solution.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
            system_rhs.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);

            hanging_node_constraints.clear ();
            DoFTools::make_hanging_node_constraints (dof_handler,
                    hanging_node_constraints);
            hanging_node_constraints.close ();
        }


    template <int dim>
        void ElasticProblem<dim>::assemble_system ()
        {
            QGauss<dim>  quadrature_formula(2);
            FEValues<dim> fe_values (fe, quadrature_formula,
                    update_values   | update_gradients |
                    update_quadrature_points | update_JxW_values);

            const unsigned int   dofs_per_cell = fe.dofs_per_cell;
            const unsigned int   n_q_points    = quadrature_formula.size();

            FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
            Vector<double>       cell_rhs (dofs_per_cell);

            std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

            std::vector<double>     lambda_values (n_q_points);
            std::vector<double>     mu_values (n_q_points);

            // ConstantFunction<dim> mu(1.);
            // cdbl E_1   = 1.0;  cdbl E_2   = 10.0;
            // // cdbl E_1   = 1.0;  cdbl E_2   = 1.0;
            // cdbl pua_1 = 0.25; cdbl pua_2 = 0.25;
            
            LambdaValue<dim> lambda(E_1, pua_1, E_2, pua_2);
            MuValue<dim> mu(E_1, pua_1, E_2, pua_2);

            vec<ATools::FourthOrderTensor> C(2);
            EPTools ::set_isotropic_elascity{yung : E_1, puasson : pua_1}(C[0]);
            EPTools ::set_isotropic_elascity{yung : E_2, puasson : pua_2}(C[1]);

            enum {x, y, z};
            vec<arr<arr<typename Nikola::SourceVector<2>::Func, 2>, 2>> U(2);
            U[0][x][x] = [&C] (const dealii::Point<2> &p) {return C[0][x][x][z][z];}; //Uz
            U[0][x][y] = [&C] (const dealii::Point<2> &p) {return 0.0;}; //Uz
            U[0][y][x] = [&C] (const dealii::Point<2> &p) {return 0.0;}; //Uz
            U[0][y][y] = [&C] (const dealii::Point<2> &p) {return C[0][x][x][z][z];}; //Uz
            U[1][x][x] = [&C] (const dealii::Point<2> &p) {return C[1][x][x][z][z];}; //Uz
            U[1][x][y] = [&C] (const dealii::Point<2> &p) {return 0.0;}; //Uz
            U[1][y][x] = [&C] (const dealii::Point<2> &p) {return 0.0;}; //Uz
            U[1][y][y] = [&C] (const dealii::Point<2> &p) {return C[1][x][x][z][z];}; //Uz

            // vec<arr<typename Nikola::SourceVector<2>::Func, 2>> tau(2);
            // tau[0][x] = [] (const dealii::Point<2> &p) {return 0.0;};
            // tau[0][y] = [] (const dealii::Point<2> &p) {return 0.0;};
            // tau[1][x] = [] (const dealii::Point<2> &p) {return 0.0;};
            // tau[1][y] = [] (const dealii::Point<2> &p) {return 0.0;};

            SourceVectorFeature<2> element_rhsv1 (U, dof_handler.get_fe());
            // SourceVectorPolyMaterials<2> element_rhsv2 (tau, dof_handler.get_fe());
            std::map<u32, dbl> list_boundary_values;

            RightHandSide<dim> right_hand_side;
            std::vector<Vector<double> > rhs_values (n_q_points,
                    Vector<double>(dim));


            typename DoFHandler<dim>::active_cell_iterator
                cell = dof_handler.begin_active(),
                     endc = dof_handler.end();
            for (; cell!=endc; ++cell)
                if (cell->subdomain_id() == this_mpi_process)
                {
                    cell_matrix = 0;
                    cell_rhs = 0;

                    if (cell->center()(1) > 0.5)
                    {
                        cell->set_material_id(1);
                        
                    };

                    fe_values.reinit (cell);

                    lambda.value_list (fe_values.get_quadrature_points(), lambda_values);
                    mu.value_list     (fe_values.get_quadrature_points(), mu_values);

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        const unsigned int
                            component_i = fe.system_to_component_index(i).first;

                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            const unsigned int
                                component_j = fe.system_to_component_index(j).first;

                            for (unsigned int q_point=0; q_point<n_q_points;
                                    ++q_point)
                            {
                                cell_matrix(i,j)
                                    +=
                                    (
                                     (fe_values.shape_grad(i,q_point)[component_i] *
                                      fe_values.shape_grad(j,q_point)[component_j] *
                                      lambda_values[q_point])
                                     +
                                     (fe_values.shape_grad(i,q_point)[component_j] *
                                      fe_values.shape_grad(j,q_point)[component_i] *
                                      mu_values[q_point])
                                     +
                                     ((component_i == component_j) ?
                                      (fe_values.shape_grad(i,q_point) *
                                       fe_values.shape_grad(j,q_point) *
                                       mu_values[q_point])  :
                                      0)
                                    )
                                    *
                                    fe_values.JxW(q_point);
                            }
                        }
                    }

                    right_hand_side.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

                    element_rhsv1 .update_on_cell (cell);
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        cell_rhs(i) = element_rhsv1(i);
                        // std::cout << cell_rhs(i) << std::endl;
                        // const unsigned int
                        //     component_i = fe.system_to_component_index(i).first;
                        //
                        // for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                        //     cell_rhs(i) += fe_values.shape_value(i,q_point) *
                        //         rhs_values[q_point](component_i) *
                        //         fe_values.JxW(q_point);
                    }

                    cell->get_dof_indices (local_dof_indices);
                    hanging_node_constraints
                        .distribute_local_to_global(cell_matrix, cell_rhs,
                                local_dof_indices,
                                system_matrix, system_rhs);

                    for (st i = 0; i < 8; ++i)
                    {
                        if (std::abs(cell->vertex(i)(x) - 5.0 ) < 1.0e-5)
                        {
                            cst v = cell->vertex_dof_index(i, x);
                            if (list_boundary_values.find(v) == list_boundary_values.end())
                            {
                                list_boundary_values.insert(std::pair<u32, dbl>(v, 0.0));
                                // list_boundary_values.emplace(v, 0.0);
                            }; 
                        };
                    };
                }

            system_matrix.compress(VectorOperation::add);
            system_rhs.compress(VectorOperation::add);
            EPTools ::print_move<2> (system_rhs, dof_handler, "rhsv.gpd");

            // std::map<types::global_dof_index,double> boundary_values;
            // VectorTools::interpolate_boundary_values (dof_handler,
            //         0,
            //         ZeroFunction<dim>(dim),
            //         boundary_values);
            // MatrixTools::apply_boundary_values (boundary_values,
            //         system_matrix, solution,
            //         system_rhs, false);
            // EPTools ::print_move<2> (system_rhs, dof_handler, "rhsv2.gpd");
            dealii::MatrixTools::apply_boundary_values (
                    list_boundary_values,
                    system_matrix,
                    solution,
                    system_rhs);
        }


    template <int dim>
        unsigned int ElasticProblem<dim>::solve ()
        {
            SolverControl solver_control (100000, 1e-8*system_rhs.l2_norm());
            PETScWrappers::SolverCG cg (solver_control, mpi_communicator);

            PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

            // std::cout << "11111" << std::endl;
            cg.solve (system_matrix, solution, system_rhs,
                    preconditioner);
            // std::cout << "22222" << std::endl;

            Vector<double> localized_solution (solution);
            // std::cout << "22223" << std::endl;

            hanging_node_constraints.distribute (localized_solution);
            // std::cout << "22224" << std::endl;

            solution = localized_solution;
            // std::cout << "22225" << std::endl;

            return solver_control.last_step();
        }

    template <int dim>
        void ElasticProblem<dim>::refine_grid ()
        {
            const Vector<double> localized_solution (solution);

            Vector<float> local_error_per_cell (triangulation.n_active_cells());
            KellyErrorEstimator<dim>::estimate (dof_handler,
                    QGauss<dim-1>(2),
                    typename FunctionMap<dim>::type(),
                    localized_solution,
                    local_error_per_cell,
                    ComponentMask(),
                    0,
                    MultithreadInfo::n_threads(),
                    this_mpi_process);

            const unsigned int n_local_cells
                = GridTools::count_cells_with_subdomain_association (triangulation,
                        this_mpi_process);
            PETScWrappers::MPI::Vector
                distributed_all_errors (mpi_communicator,
                        triangulation.n_active_cells(),
                        n_local_cells);

            for (unsigned int i=0; i<local_error_per_cell.size(); ++i)
                if (local_error_per_cell(i) != 0)
                    distributed_all_errors(i) = local_error_per_cell(i);
            distributed_all_errors.compress (VectorOperation::insert);


            const Vector<float> localized_all_errors (distributed_all_errors);

            GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                    localized_all_errors,
                    0.3, 0.03);
            triangulation.execute_coarsening_and_refinement ();
        }


    template <int dim>
        void ElasticProblem<dim>::output_results (const unsigned int cycle) const
        {
            const Vector<double> localized_solution (solution);
            if (this_mpi_process == 0)
            {
                std::ostringstream filename;
                // filename << "solution-" << cycle << ".gmv";
                filename << "solution.gpd";

                std::ofstream output (filename.str().c_str());

                DataOut<dim> data_out;
                data_out.attach_dof_handler (dof_handler);

                std::vector<std::string> solution_names;
                switch (dim)
                {
                    case 1:
                        solution_names.push_back ("displacement");
                        break;
                    case 2:
                        solution_names.push_back ("x_displacement");
                        solution_names.push_back ("y_displacement");
                        break;
                    case 3:
                        solution_names.push_back ("x_displacement");
                        solution_names.push_back ("y_displacement");
                        solution_names.push_back ("z_displacement");
                        break;
                    default:
                        Assert (false, ExcInternalError());
                }

                data_out.add_data_vector (localized_solution, solution_names);

                std::vector<unsigned int> partition_int (triangulation.n_active_cells());
                GridTools::get_subdomain_association (triangulation, partition_int);

                const Vector<double> partitioning(partition_int.begin(),
                        partition_int.end());

                data_out.add_data_vector (partitioning, "partitioning");

                data_out.build_patches ();
                data_out.write_gnuplot (output);
                {
                    LambdaValue<dim> lambda(E_1, pua_1, E_2, pua_2);
                    MuValue<dim> mu(E_1, pua_1, E_2, pua_2);

                    vec<ATools::FourthOrderTensor> C(2);
                    EPTools ::set_isotropic_elascity{yung : E_1, puasson : pua_1}(C[0]);
                    EPTools ::set_isotropic_elascity{yung : E_2, puasson : pua_2}(C[1]);

                    arr<dealii::Vector<dbl>, 2> stress;
                    get_nikola_stress (localized_solution, dof_handler, C, stress);
                    EPTools ::print_move<2> (stress[0], dof_handler, "stress_x_1.gpd");
                    EPTools ::print_move<2> (stress[1], dof_handler, "stress_y_1.gpd");
                };
            }
        }


    template <int dim>
        void ElasticProblem<dim>::run ()
        {
            std::cout << n_mpi_processes << std::endl;

            // GridGenerator::hyper_cube (triangulation, -1, 1);
            // triangulation.refine_global (5);
            set_laminat(triangulation, [](cdbl x){return 0.5;}, 10.0, 1.0, 20);

            triangulation.n_active_cells();

            setup_system ();

            assemble_system ();
            const unsigned int n_iterations = solve ();
            std::cout << "22226" << std::endl;

            output_results (0);
        }
}


int main (int argc, char **argv)
{
    try
    {
        using namespace dealii;
        using namespace Step17;

        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
        {
            deallog.depth_console (0);

            ElasticProblem<2> elastic_problem;
            elastic_problem.run ();
        }
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
        std::cerr << "Exception on processing: " << std::endl
            << exc.what() << std::endl
            << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;

        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
        std::cerr << "Unknown exception!" << std::endl
            << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
        return 1;
    }

    return 0;
}
