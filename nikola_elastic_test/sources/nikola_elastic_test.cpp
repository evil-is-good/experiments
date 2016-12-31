#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>



#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_relaxation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/relaxation_block.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/vector_tools.h>



#include "../../feminist/calculation_core/src/blocks/general/domain/domain.h"
#include "../../feminist/calculation_core/src/blocks/general/laplacian/vector/laplacian_vector.h"
#include "../../feminist/calculation_core/src/blocks/general/source/scalar/source_scalar.h"
#include "../../feminist/calculation_core/src/blocks/general/boundary_value/boundary_value.h"
#include "../../feminist/calculation_core/src/blocks/general/assembler/assembler.h"
#include "../../feminist/calculation_core/src/blocks/general/system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../../feminist/calculation_core/src/blocks/general/additional_tools/trivial_prepare_system_equations/trivial_prepare_system_equations.h"
#include "../../feminist/calculation_core/src/blocks/general/additional_tools/apply_boundary_value/scalar/apply_boundary_value_scalar.h"
#include "../../feminist/calculation_core/src/blocks/general/geometric_tools/geometric_tools.h"


#include "../../feminist/calculation_core/src/blocks/special/elastic_problem_tools/elastic_problem_tools.h"

#include "../../feminist/calculation_core/src/blocks/general/grid_generator/grid_generator.h"
#include "../../bloks/src/grid_generator/grid_generator.h"

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

    // vec<prmt::Point<2>> border;
    // vec<st> type_border;
    // GridGenerator::give_rectangle_with_border_condition(
    //         border,
    //         type_border,
    //         arr<st, 4>{1,3,2,4},
    //         10,
    //         prmt::Point<2>(0.0, 0.0), prmt::Point<2>(lx, ly));
    // vec<vec<prmt::Point<2>>> inclusion(1);
    // for (st i = 0; i < np; ++i)
    // {
    //     cdbl X = i * lx / np;
    //     inclusion[0].push_back(prmt::Point<2>(X, delimiter(X)));
    // };
    // CGALGridGenerator::set_grid_with_constraints(triangulation, border, inclusion, type_border);

    // dealii::GridGenerator::hyper_rectangle(triangulation, 
    //         dealii::Point<2>(0.0, 0.0),
    //         dealii::Point<2>(lx, ly));

    {
        std::vector< dealii::Point< 2 > > v (8);

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
        // c .push_back (dealii::CellData<2>{{0, 1, 3, 2}, 0});
        // c .push_back (dealii::CellData<2>{{3, 5, 4, 2}, 1});

        // dealii::SubCellData b;
        // {
        //     dealii::CellData<1> tmp;
        //     tmp.vertices[0]=4;tmp.vertices[1]=2;tmp.boundary_id=0;
        //     b.boundary_lines .push_back (tmp);
        //     tmp.vertices[0]=2;tmp.vertices[1]=0;tmp.boundary_id=0;
        //     b.boundary_lines .push_back (tmp);
        //     tmp.vertices[0]=0;tmp.vertices[1]=1;tmp.boundary_id=1;
        //     b.boundary_lines .push_back (tmp);
        //     tmp.vertices[0]=1;tmp.vertices[1]=3;tmp.boundary_id=2;
        //     b.boundary_lines .push_back (tmp);
        //     tmp.vertices[0]=3;tmp.vertices[1]=5;tmp.boundary_id=2;
        //     b.boundary_lines .push_back (tmp);
        //     tmp.vertices[0]=5;tmp.vertices[1]=4;tmp.boundary_id=3;
        //     b.boundary_lines .push_back (tmp);
        // };
        // b.boundary_lines .push_back (dealii::CellData<1>{4, 2, 0});
        // b.boundary_lines .push_back (dealii::CellData<1>{2, 0, 0});
        // b.boundary_lines .push_back (dealii::CellData<1>{0, 1, 1});
        // b.boundary_lines .push_back (dealii::CellData<1>{1, 3, 2});
        // b.boundary_lines .push_back (dealii::CellData<1>{3, 5, 2});
        // b.boundary_lines .push_back (dealii::CellData<1>{5, 4, 3});

        // dealii::GridReordering<2> ::reorder_cells (c);
        triangulation .create_triangulation (v, c, dealii::SubCellData());
    };
    triangulation.refine_global(3);
    // for (st i = 0; i < 3; ++i)
    // {
    //     {
    //         auto cell = triangulation .begin_active();
    //         auto end_cell = triangulation .end();
    //         for (; cell != end_cell; ++cell)
    //         {
    //             auto c = cell->center();
    //             if ((c(x) < 1.0) or (c(x) > 9.0))
    //             {
    //                 // cell->set_refinement_case (
    //                 //         dealii::RefinementCase<2>(dealii::RefinementPossibilities<2>::Possibilities::cut_x));
    //                 dealii::RefinementCase<2> ref_case(1);
    //                 cell->set_refine_flag (ref_case);
    //             };
    //         };
    //     };
    //     triangulation.execute_coarsening_and_refinement();
    // };

    std::ofstream out ("grid-igor.eps");
    dealii::GridOut grid_out;
    grid_out.write_eps (triangulation, out);

    {
        auto cell = triangulation .begin_active();
        auto end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            auto c = cell->center();
            // dealii::Point<2> midle_p(0.0, 0.0);
            //
            // for (size_t i = 0; i < 4; ++i)
            // {
            //     midle_p(0) += cell->vertex(i)(0);
            //     midle_p(1) += cell->vertex(i)(1);
            // };
            // midle_p(0) /= 4.0;
            // midle_p(1) /= 4.0;

            // printf("%f %f\n", midle_p(0), midle_p(1));

            if (c(y) > delimiter(c(x)))
            {
                // std::cout << c(y) << " " <<  delimiter(c(x)) << " " <<  1 << std::endl;
                cell->set_material_id(1);
            }
            else
            {
                // std::cout << c(y) << " " <<  delimiter(c(x)) << " " <<  0 << std::endl;
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

void get_nikola_mean_stress (const dealii::Vector<dbl> &move, 
        const dealii::DoFHandler<2> &dof_handler,
        const vec<ATools::FourthOrderTensor> &C,
        arr<dealii::Vector<dbl>, 2> &stress)
{
    enum {x, y, z};

    cu32 dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    dealii::QGauss<2> quadrature_formula(2);

    dealii::FEValues<2> fe_values (dof_handler.get_fe(), quadrature_formula,
            dealii::update_gradients | dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

    cst num_quad_points = quadrature_formula.size();

    stress[0] .reinit (dof_handler.n_dofs());
    stress[1] .reinit (dof_handler.n_dofs());

    vec<st> N(dof_handler.n_dofs());

    for (
            auto cell = dof_handler.begin_active(); 
            cell     != dof_handler.end(); 
            ++cell
        )
    {
        fe_values .reinit (cell);

        dbl cell_area = 0.0;
        for (st q_point = 0; q_point < num_quad_points; ++q_point)
        {
            cell_area += fe_values.JxW(q_point);
        };

        arr<arr<dbl, 2>, 2> grad;
        for (st i = 0; i < 2; ++i)
        {
            for (st j = 0; j < 2; ++j)
            {
                grad[i][j] = 0.0;
                // std::cout << grad[i][j] << std::endl;
            };
        };
        for (st point = 0; point < 4; ++point)
        {
            for (st q_point = 0; q_point < num_quad_points; ++q_point)
            {
                for (st i = 0; i < 2; ++i)
                {
                    for (st j = 0; j < 2; ++j)
                    {
                        grad[i][j] += 
                            move(cell->vertex_dof_index (point, i)) *
                            fe_values.shape_grad (point*2+i, q_point)[j] *
                            fe_values.JxW(q_point);
                        // std::cout <<
                        //     move(cell->vertex_dof_index (point, i)) << " " <<  
                        //     fe_values.shape_grad (point*2+i, q_point)[j] << " " <<  
                        //     fe_values.JxW(q_point) << " " <<
                        //     move(cell->vertex_dof_index (point, i)) *
                        //     fe_values.shape_grad (point*2+i, q_point)[j] *
                        //     fe_values.JxW(q_point) << " " <<  
                        //     grad[i][j]
                        //     << std::endl;
                    };
                };
            };
        };
        for (st i = 0; i < 2; ++i)
        {
            for (st j = 0; j < 2; ++j)
            {
                grad[i][j] /= cell_area;// * 4.0;
                // std::cout << grad[i][j] << std::endl;
            };
        };
        for (st i = 0; i < 4; ++i)
        {
            auto mat_id = cell->material_id();
            auto indx_1 = cell->vertex_dof_index(i, x);
            auto indx_2 = cell->vertex_dof_index(i, y);

            stress[x][indx_1] += 
                // grad[x][x];//C[mat_id][0][0][0][0];
                C[mat_id][x][x][x][x] * grad[x][x] +
                C[mat_id][x][x][x][y] * grad[x][y] +
                C[mat_id][x][x][y][x] * grad[y][x] +
                C[mat_id][x][x][y][y] * grad[y][y] + C[mat_id][x][x][z][z];
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
        // break;
    };
    for (st i = 0; i < N.size(); ++i)
    {
        stress[0][i] /= N[i];
        stress[1][i] /= N[i];
    };
};

void solve_nikola_elastic_problem (dbl E_1, cdbl pua_1, cdbl E_2, cdbl pua_2, Domain<2>& domain)
{
    enum {x, y, z};
    dealii::FESystem<2,2> fe 
        (dealii::FE_Q<2,2>(1), 2);
    domain.dof_init (fe);

    SystemsLinearAlgebraicEquations slae;
    ATools ::trivial_prepare_system_equations (slae, domain);

    LaplacianVector<2> element_matrix (domain.dof_handler.get_fe());
    element_matrix.C .resize (2);
    EPTools ::set_isotropic_elascity{yung : E_1, puasson : pua_1}(element_matrix.C[0]);
    EPTools ::set_isotropic_elascity{yung : E_2, puasson : pua_2}(element_matrix.C[1]);

    // T2.2
    vec<arr<arr<typename Nikola::SourceVector<2>::Func, 2>, 2>> U(2);
    U[0][x][x] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z];}; //Uz
    U[0][x][y] = [&element_matrix] (const dealii::Point<2> &p) {return 0.0;}; //Uz
    U[0][y][x] = [&element_matrix] (const dealii::Point<2> &p) {return 0.0;}; //Uz
    U[0][y][y] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z];}; //Uz
    U[1][x][x] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z];}; //Uz
    U[1][x][y] = [&element_matrix] (const dealii::Point<2> &p) {return 0.0;}; //Uz
    U[1][y][x] = [&element_matrix] (const dealii::Point<2> &p) {return 0.0;}; //Uz
    U[1][y][y] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z];}; //Uz

    vec<arr<typename Nikola::SourceVector<2>::Func, 2>> tau(2);
    tau[0][x] = [] (const dealii::Point<2> &p) {return 0.0;};
    tau[0][y] = [] (const dealii::Point<2> &p) {return 0.0;};
    tau[1][x] = [] (const dealii::Point<2> &p) {return 0.0;};
    tau[1][y] = [] (const dealii::Point<2> &p) {return 0.0;};

    SourceVectorFeature<2> element_rhsv1 (U, domain.dof_handler.get_fe());
    SourceVectorPolyMaterials<2> element_rhsv2 (tau, domain.dof_handler.get_fe());

    Assembler ::assemble_matrix<2> (slae.matrix, element_matrix, domain.dof_handler);

    Assembler ::assemble_rhsv<2> (slae.rhsv, element_rhsv1, domain.dof_handler);
    Assembler ::assemble_rhsv<2> (slae.rhsv, element_rhsv2, domain.dof_handler);
    EPTools ::print_move<2> (slae.rhsv, domain.dof_handler, "rhsv.gpd");

    {
        std::map<u32, dbl> list_boundary_values;
        auto cell = domain.dof_handler .begin_active();
        auto end_cell = domain.dof_handler .end();
        for (; cell != end_cell; ++cell)
        {
            for (st i = 0; i < 8; ++i)
            {
                // for (st ort = 0; ort < 3; ++ort)
                // {
                    if (std::abs(cell->vertex(i)(x) - 5.0 ) < 1.0e-5)
                    {
                        cst v = cell->vertex_dof_index(i, x);
                        if (list_boundary_values.find(v) == list_boundary_values.end())
                        {
                            // list_boundary_values.insert(std::pair<u32, dbl>(v, 0.0));
                            list_boundary_values.emplace(v, 0.0);
                        }; 
                    };
                    // if ((std::abs(cell->vertex(i)(x) - 0.5) < 1.0e-5) and (std::abs(cell->vertex(i)(y)) < 1.0e-5))
                    // {
                    //     cst v = cell->vertex_dof_index(i, y);
                    //     if (list_boundary_values.find(v) == list_boundary_values.end())
                    //     {
                    //         list_boundary_values.insert(std::pair<u32, dbl>(v, 0.0));
                    //     };
                    // };

                    // if (std::abs(cell->vertex(i)(x) - 5.0 ) < 1.0e-5)
                    // {
                    //     cst v = cell->vertex_dof_index(i, y);
                    //     if (list_boundary_values.find(v) == list_boundary_values.end())
                    //     {
                    //         list_boundary_values.insert(std::pair<u32, dbl>(v, 0.0));
                    //     }; 
                    // };
                // };
            };
        };

        dealii::MatrixTools::apply_boundary_values (
                list_boundary_values,
                slae.matrix,
                slae.solution,
                slae.rhsv);
    };

    dealii::SolverControl solver_control (1000000, 1e-8);
    dealii::SolverCG<> solver (solver_control);
    solver.solve (
            slae.matrix,
            slae.solution,
            slae.rhsv
            ,dealii::PreconditionIdentity()
            );

    EPTools ::print_move<2> (slae.solution, domain.dof_handler, "move.gpd");
    {
        arr<dealii::Vector<dbl>, 2> stress;
        get_nikola_stress (slae.solution, domain.dof_handler, element_matrix.C, stress);
        EPTools ::print_move<2> (stress[x], domain.dof_handler, "stress_x_1.gpd");
        EPTools ::print_move<2> (stress[y], domain.dof_handler, "stress_y_1.gpd");
    };
    {
        arr<dealii::Vector<dbl>, 2> stress;
        get_nikola_mean_stress (slae.solution, domain.dof_handler, element_matrix.C, stress);
        EPTools ::print_move<2> (stress[x], domain.dof_handler, "stress_x_2.gpd");
        EPTools ::print_move<2> (stress[y], domain.dof_handler, "stress_y_2.gpd");
    {
        dbl max = 0.0;
        for (st i = 0; i < stress[y].size(); ++i)
        {
            if (std::abs(max) < std::abs(stress[y](i)))
            {
                max = stress[y](i);
            };
        };
        std::cout << "\x1B[36m max = " << max << "\x1B[0m     File: " << __FILE__ << " Line: " << __LINE__ << std::endl; //DEBAG OUT
    };
    };
};

int main()
{
    enum {x, y, z};
    Domain<2> domain;
    // dealii::GridIn<2> gridin;
    // gridin.attach_triangulation(domain.grid);
    // // std::ifstream f("../sources/tst1.msh");
    // std::ifstream f("tst1_true.msh");
    // // // std::ifstream f("../sources/circle_R2_R1.msh");
    // gridin.read_msh(f);
    // std::ofstream out ("grid-igor.eps");
    // dealii::GridOut grid_out;
    // grid_out.write_eps (domain.grid, out);
    set_laminat(domain.grid, [](cdbl x){return 0.5;}, 10.0, 1.0, 20);
    // set_laminat(tria, [](cdbl x){return std::sin(x*M_PI)/4.0+0.5;}, 2.0, 1.0, 20);
    // dealii::GridGenerator::hyper_cube(domain.grid, 0.0, 1.0);
    // domain.grid.refine_global(3);
    // solve_nikola_elastic_problem(1.0, 0.25, 10.0, 0.25,  domain);
    // solve_nikola_elastic_problem(1.0, 0.2, 10.0, 0.1,  domain);
    solve_nikola_elastic_problem(1.0, 0.2, 1.0, 0.1,  domain);
    // solve_nikola_elastic_problem(1.0, 0.1, 1.0, 0.1, domain);
    return 0;
}
