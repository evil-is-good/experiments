#include "pch.h"

#include "../../feminist/calculation_core/src/blocks/general/grid_generator/grid_generator.h"
#include "../../bloks/src/grid_generator/grid_generator.h"

#include "../../feminist/calculation_core/src/blocks/special/elastic_problem_tools/elastic_problem_tools.h"

// #include "../../../calculation_core/src/blocks/general/additional_tools/apply_boundary_value/vector/apply_boundary_value_vector.h"
// #include "../../../calculation_core/src/blocks/general/grid_generator/grid_generator.h"
//
// #include "../../../calculation_core/src/blocks/special/heat_conduction_problem_tools/heat_conduction_problem_tools.h"
// #include "../../../calculation_core/src/blocks/special/elastic_problem_tools/elastic_problem_tools.h"

// void set_ball_true_ordered_speciment(
//         dealii::Triangulation<3> &triangulation,
//         cdbl R,
//         cst n_ref_cells,
//         cst n_refine,
//         cst n_slices,
//         cdbl H)
// {
//     enum {x, y, z};
//     cdbl eps = 1.0e-8;
//
//     cst n_cells = 1 << n_ref_cells;
//     cst all_n_cells = n_cells * n_cells * n_slices;
//     cdbl start_center = 1.0 / (2 << n_ref_cells);
//     cdbl start_center_z = (1.0 / n_slices) / 2.0;
//     cdbl size_cell = 1.0 / n_cells;
//     cdbl size_cell_z = 1.0 / n_slices;
//     cdbl Rc = R / n_cells;
//     vec<dealii::Point<2>> center(all_n_cells);
//     {
//         st count = 0;
//         dealii::Point<3> c(start_center, start_center, start_center_z);
//         for (st i = 0; i < n_cells; ++i)
//         {
//             c(y) = start_center;
//             for (st j = 0; j < n_cells; ++j)
//             {
//                 c(z) = start_center_z;
//                 for (st k = 0; k < n_slices; ++k)
//                 {
//                     center[count] = c;
//                     c(z) += size_celli_z;
//                     ++count;
//                 };
//                 c(y) += size_cell;
//             };
//             c(x) += size_cell;
//         };
//     };
//     std::cout << "start_center = " << start_center << std::endl; //DEBAG OUT
//     std::cout << "size_cell = " << size_cell << std::endl; //DEBAG OUT
//     std::cout << "Rc = " << Rc << std::endl; //DEBAG OUT    
//     std::cout << "n_cells = " << n_cells << std::endl; //DEBAG OUT
//
//     dealii::GridGenerator ::hyper_rectangle (triangulation, 
//             dealii::Point<3>(0.0, 0.0, 0.0), 
//             dealii::Point<3>(1.0, 1.0, H));
//     triangulation .refine_global (n_refine);
//
//     cdbl fe_size = 1.0 / (1 << n_refine);
//     cdbl max_dist = std::sqrt(std::pow(fe_size, 2.0) * 3.0) / 2.0;
//
//     {
//         auto cell = triangulation .begin_active();
//         auto end_cell = triangulation .end();
//         for (; cell != end_cell; ++cell)
//         {
//             for (st i = 0; i < dealii::GeometryInfo<3>::vertices_per_cell; ++i)
//             {
//                 auto p = cell->vertex(i);
//                 cdbl r = center.distance(p); 
//                 cdbl dist = std::abs(r - radius);
//                 // cdbl dist = std::abs(center.distance(p) - R); 
//                 if (dist < max_dist)
//                 {
//                     cell->vertex(i) -= center;
//                     cell->vertex(i) *= radius/r;
//                     cell->vertex(i) += center;
//                 };
//             };
//         };
//     };
//
//     {
//         auto cell = triangulation .begin_active();
//         auto end_cell = triangulation .end();
//         for (; cell != end_cell; ++cell)
//         {
//             if (center.distance(cell->center()) < radius)
//             {
//                 cell->set_material_id(1);
//             }
//             else
//                 cell->set_material_id(0);
//         };
//     };
//
//     
//
//
//     dealii::Triangulation<2> tria2d;
//     // dealii::GridGenerator ::hyper_cube (tria2d, -0.5, 0.5);
//     dealii::GridGenerator ::hyper_cube (tria2d, 0.0, 1.0);
//     tria2d .refine_global (level_density);
//     // tria2d .refine_global (1);
//     std::cout << "n_cells " << tria2d.n_active_cells() << std::endl;
//
//     cdbl max_dist = std::sqrt(std::pow(1.0 / std::pow(2.0, level_density), 2.0) * 2.0) / 2.0;
//
//     auto face = tria2d.begin_active_face();
//     auto endf = tria2d.end_face();
//     for (; face != endf; ++face)
//     {
//         for (st i = 0; i < 2; ++i)
//         {
//             auto p = face->vertex(i);
//             for (st j = 0; j < all_n_cells; ++j)
//             {
//                 cdbl r = center[j].distance(p); 
//                 cdbl dist = std::abs(r - Rc);
//                 if (dist < max_dist)
//                 {
//                     face->vertex(i) -= center[j];
//                     face->vertex(i) *= Rc/r;
//                     face->vertex(i) += center[j];
//                 };
//             };
//         };
//     };
//
//     {
//         auto cell = tria2d .begin_active();
//         auto end_cell = tria2d .end();
//         for (; cell != end_cell; ++cell)
//         {
//             cell->set_material_id(0);
//             for (st i = 0; i < all_n_cells; ++i)
//             {
//                 if (center[i].distance(cell->center()) < Rc)
//                     // if (
//                     //         (std::abs(center[i](x) - cell->center()(x)) < Rc) and
//                     //         (std::abs(center[i](y) - cell->center()(y)) < Rc)
//                     //         )
//                 {
//                     cell->set_material_id(1);
//                 };
//             };
//         };
//     };
//
//     {
//         std::ofstream out ("grid-speciment.eps");
//         dealii::GridOut grid_out;
//         grid_out.write_eps (tria2d, out);
//     };
//
//     extrude_triangulation (tria2d, n_slices, H, triangulation);
//
//     // {
//     //     dealii::Triangulation<3>::active_cell_iterator
//     //         cell = triangulation .begin_active(),
//     //              end_cell = triangulation .end();
//     //     for (; cell != end_cell; ++cell)
//     //     {
//     //         dealii::Point<2> midle_p(0.0, 0.0);
//     //
//     //         for (size_t i = 0; i < 8; ++i)
//     //         {
//     //             st count = 0;
//     //             for (st j = 0; j < 2; ++j)
//     //             {
//     //                 midle_p(count) += cell->vertex(i)(j);
//     //             };
//     //         };
//     //         midle_p(0) /= 8.0;
//     //         midle_p(1) /= 8.0;
//     //
//     //         cell->set_material_id(1);
//     //         for (st i = 0; i < all_n_cells; ++i)
//     //         {
//     //             if (center[i].distance(midle_p) < Rc)
//     //             {
//     //                 cell->set_material_id(0);
//     //             };
//     //         };
//     //     };
//     // };
// 

void gen_laminate(dealii::Triangulation< 3 > &triangulation, 
        cdbl W, cst n_slices, cst n_refine)
{
    cdbl width = W / n_slices;
    std::cout << width << std::endl;
    dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
    triangulation .refine_global (n_refine);
    {
        dealii::Triangulation<3>::active_cell_iterator
            cell = triangulation .begin_active(),
                 end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            for (st i = 0; i < n_slices; ++i)
            {
                if ((i*2*width < cell->center()(0)) and (cell->center()(0) < (i*2+1)*width))
                {
                    // std::cout << cell->center()(0) << std::endl;
                    cell->set_material_id(1);
                    break;
                }
                else
                if (((i*2+1)*width < cell->center()(0)) and (cell->center()(0) < (i*2+2)*width))
                {
                    cell->set_material_id(0);
                    break;
                };
            };
        };
    };
};

void gen_laminate(dealii::Triangulation< 3 > &triangulation, 
        cst n_slices, cst n_refine, cdbl H, cst n_v_slices)
{
    cdbl width = 1.0 / n_slices;
    arr<cdbl, 4> mark = {0.0, (1.0/4.0)/n_slices, (3.0/4.0)/n_slices, 1.0};
    // dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
    // triangulation .refine_global (n_refine);
    dealii::Triangulation<2> tria2d;
    dealii::GridGenerator ::hyper_cube (tria2d, 0.0, 1.0);
    tria2d .refine_global (n_refine);
    {
        auto cell = tria2d .begin_active();
        auto end_cell = tria2d .end();
        for (; cell != end_cell; ++cell)
        {
            for (st i = 0; i < n_slices; ++i)
            {
                cdbl offset = i * width;
                if ((mark[0] + offset < cell->center()(0)) and (cell->center()(0) < mark[1] + offset))
                {
                    cell->set_material_id(0);
                };
                if ((mark[1] + offset < cell->center()(0)) and (cell->center()(0) < mark[2] + offset))
                {
                    cell->set_material_id(1);
                };
                if ((mark[2] + offset < cell->center()(0)) and (cell->center()(0) < mark[3] + offset))
                {
                    cell->set_material_id(0);
                };
            };
        };
    };
    GridGenerator::extrude_triangulation (tria2d, n_v_slices, H, triangulation);
};

arr<dbl, 3> calc_sum_move(
        dealii::DoFHandler<3>& dof_h, 
        dealii::Vector<dbl>& solution, 
        const arr<dbl, 3>& size)
{
    enum {x, y, z};
    cdbl eps = 1.0e-8;
    arr<dbl, 3> sum = {0.0, 0.0, 0.0};

    {
        dealii::QGauss<2>  quadrature_formula(2);

        dealii::FEFaceValues<3> fe_values (dof_h.get_fe(), quadrature_formula,
                dealii::update_values | dealii::update_quadrature_points | dealii::update_JxW_values);

        cst n_q_points = quadrature_formula.size();
        cst dofs_per_cell = 8;//8;//fe_values.get_fe().dofs_per_cell;

        auto cell = dof_h.begin_active();
        auto end_cell = dof_h.end();
        for (; cell != end_cell; ++cell)
        {
            if (cell->at_boundary())
            {
                for (st i = 0; i < 3; ++i)
                {
                    for (st j = 0; j < 6; ++j)
                    {
                        if (std::abs(cell->face(j)->center()(i) - size[i]) < eps)
                        {
                            fe_values .reinit (cell, j);

                            dbl integ = 0.0;
                            for (st q_point = 0; q_point < n_q_points; ++q_point)
                            {
                                dbl f_q = 0.0;
                                for (st n = 0; n < dofs_per_cell; ++n)
                                {
                                    cst dof = 3 * n + i;
                                    f_q += solution(cell->vertex_dof_index(n, i)) * fe_values.shape_value(dof, q_point);
                                };
                                integ += f_q * fe_values.JxW(q_point);
                            };
                            sum[i] += integ;
                            break;
                        };
                    };
                };
            };
        };
    };

    std::cout << sum[x] << " " <<  sum[y] << " " <<  sum[z] << std::endl;
    sum[x] /= size[y] * size[z];
    sum[y] /= size[x] * size[z];
    sum[z] /= size[x] * size[y];
    std::cout << sum[x] << " " <<  sum[y] << " " <<  sum[z] << std::endl;

    return sum;
};

arr<dbl, 3> solve_cpeciment_extension_test (
        cst direction, 
        const arr<dbl, 3>& size, 
        cdbl P,
        vec<ATools::FourthOrderTensor>& C, 
        Domain<3>& domain)
{
    enum {x, y, z};

    cst d = direction;

    dealii::FESystem<3,3> fe(dealii::FE_Q<3,3>(1), 3);
    domain.dof_init (fe);

    SystemsLinearAlgebraicEquations slae;
    ATools ::trivial_prepare_system_equations (slae, domain);

    LaplacianVector<3> element_matrix (domain.dof_handler.get_fe());
    element_matrix.C .resize (C.size());

    for (st i = 0; i < C.size(); ++i)
    {
        ATools:: copy_tensor (C[i], element_matrix.C[i]);
    };

    arr<std::function<dbl (const dealii::Point<3>&)>, 3> func {
        [] (const dealii::Point<3>) {return 0.0;},
            [] (const dealii::Point<3>) {return 0.0;},
            [] (const dealii::Point<3>) {return 0.0;}
    };
    SourceVector<3> element_rhsv (func, domain.dof_handler.get_fe());

    Assembler ::assemble_matrix<3> (slae.matrix, element_matrix, domain.dof_handler);

    // прикладываем расстягивающие усилия
    {
        auto face = domain.grid .begin_active_face();
        auto end_face = domain.grid .end_face();
        for (; face != end_face; ++face)
        {
            if (std::abs(face->center()(d)) < 1.0e-5)
            {
                face->set_boundary_id(2);
            };
            if (std::abs(face->center()(d) - size[d]) < 1.0e-5)
            {
                face->set_boundary_id(1);
            };
        };
    };

    vec<BoundaryValueVector<3>> bound (2);
    bound[0].function      = [P, d] (const dealii::Point<3> &p) {
        arr<dbl, 3> value = {0.0, 0.0, 0.0};
        value[d] = P;
        return value;};
    bound[0].boundary_id   = 1;
    bound[0].boundary_type = TBV::Neumann;
    // bound[0].boundary_type = TBV::Dirichlet;
    bound[1].function      = [P, d] (const dealii::Point<3> &p) {
        arr<dbl, 3> value = {0.0, 0.0, 0.0};
        value[d] = -P;
        return value;};
    bound[1].boundary_id   = 2;
    bound[1].boundary_type = TBV::Neumann;
    // bound[1].boundary_type = TBV::Dirichlet;

    for (auto b : bound)
        ATools ::apply_boundary_value_vector<3> (b) .to_slae (slae, domain);

    //Зажимаем точки в середине по оси z и по ося x y
    {
        std::map<u32, dbl> list_boundary_values;
        auto cell = domain.dof_handler .begin_active();
        auto end_cell = domain.dof_handler .end();
        for (; cell != end_cell; ++cell)
        {
            for (st i = 0; i < 8; ++i)
            {
                for (st ort = 0; ort < 3; ++ort)
                {
                    if (std::abs(cell->vertex(i)(ort) - size[ort] / 2.0 ) < 1.0e-5)
                    {
                        cst v = cell->vertex_dof_index(i, ort);
                        if (list_boundary_values.find(v) == list_boundary_values.end())
                        {
                            list_boundary_values.insert(std::pair<u32, dbl>(v, 0.0));
                        }; 
                    };
                };
            };
        };

        dealii::MatrixTools::apply_boundary_values (
                list_boundary_values,
                slae.matrix,
                slae.solution,
                slae.rhsv);
    };

    dealii::SolverControl solver_control (1000000, 1e-12);
    dealii::SolverCG<> solver (solver_control);
    solver.solve (
            slae.matrix,
            slae.solution,
            slae.rhsv
            ,dealii::PreconditionIdentity()
            );

    arr<str, 3> ort = {"x", "y", "z"};
    EPTools ::print_move<3> (slae.solution, domain.dof_handler, ort[d] + "/move.gpd");
    EPTools ::print_move<3> (slae.solution, domain.dof_handler, ort[d] + "/move.vtk", dealii::DataOutBase::vtk);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_x_05.gpd", x, size[x] / 2.0);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_x_1.gpd", x, size[x]);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_y_05.gpd", y, size[y] / 2.0);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_y_1.gpd", y, size[y]);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_z_05.gpd", z, size[z] / 2.0);
    EPTools ::print_move_slice (slae.solution, domain.dof_handler, ort[d] + "/move_slice_z_1.gpd", z, size[z]);

    return calc_sum_move(domain.dof_handler, slae.solution, size);
    // return arr<dbl, 3>{0.0, 0.0, 0.0};
};

    void gen_circle_in_square_true_ordered(
            dealii::Triangulation<2> &tria2d,
            cdbl R,
            cst level_density)
    {
        enum {x, y, z};
        cdbl eps = 1.0e-8;

        dealii::Point<2> center(0.5, 0.5);

        dealii::GridGenerator ::hyper_cube (tria2d, 0.0, 1.0);
        tria2d .refine_global (level_density);

        cdbl max_dist = std::sqrt(std::pow(1.0 / std::pow(2.0, level_density), 2.0) * 2.0) / 2.0;

        auto face = tria2d.begin_active_face();
        auto endf = tria2d.end_face();
        for (; face != endf; ++face)
        {
            for (st i = 0; i < 2; ++i)
            {
                auto p = face->vertex(i);
                cdbl r = center.distance(p); 
                cdbl dist = std::abs(r - R);
                if (dist < max_dist)
                {
                    face->vertex(i) -= center;
                    face->vertex(i) *= R/r;
                    face->vertex(i) += center;
                };
            };
        };

        {
            auto cell = tria2d .begin_active();
            auto end_cell = tria2d .end();
            for (; cell != end_cell; ++cell)
            {
                cell->set_material_id(0);
                if (center.distance(cell->center()) < R)
                    // if (
                    //         (std::abs(center[i](x) - cell->center()(x)) < Rc) and
                    //         (std::abs(center[i](y) - cell->center()(y)) < Rc)
                    //         )
                {
                    cell->set_material_id(1);
                };
            };
        };

        {
            std::ofstream out ("grid.eps");
            dealii::GridOut grid_out;
            grid_out.write_eps (tria2d, out);
        };
    };

int main()
{
    enum {x, y, z};

    cdbl E1 = 0.6;
    cdbl E2 = 60.0;
    // cdbl E1 = 1.0;
    // cdbl E2 = 1.0;
    vec<ATools::FourthOrderTensor> C(2);
    EPTools ::set_isotropic_elascity{yung : E1, puasson : 0.35}(C[0]);
    EPTools ::set_isotropic_elascity{yung : E2, puasson : 0.20}(C[1]);
    
    Domain<3> domain;
    cst n_slices = 5;
    cst n_ref = 4;
    cst n_p = 16;
    cst n_cell_ref = 2;
    cst n_cells = 32;
    cdbl R = sqrt(0.7 / M_PI) / 4.0;//0.45;
    std::cout << "\x1B[36m R = " << R << "\x1B[0m     File: " << __FILE__ << " Line: " << __LINE__ << std::endl; //DEBAG OUT
    cdbl H = 10.0;
    cdbl W = n_cells/4.0;
    // GridGenerator::gen_cylinder_true_ordered_speciment(domain.grid, R, n_cell_ref, n_ref, n_slices, H);
    vec<dealii::Point<2>> center;
    for (st i = 0; i < n_cells; ++i)
    {
        // center .push_back(dealii::Point<2>(0.25+i*0.5, 0.25));
        // center .push_back(dealii::Point<2>(0.25+i*0.5, 0.75));
        center .push_back(dealii::Point<2>(0.125+i*0.25, 0.125));
        center .push_back(dealii::Point<2>(0.125+i*0.25, 0.375));
        center .push_back(dealii::Point<2>(0.125+i*0.25, 0.625));
        center .push_back(dealii::Point<2>(0.125+i*0.25, 0.875));
    };

    const arr<dbl, 3> size = {W, 1.0, H};
    GridGenerator::set_cylinder_in_rectangular_cgal(domain.grid, size, center, R, n_p, n_slices);
    // GridGenerator::gen_cylinder_true_ordered_speciment_for_x_extension(domain.grid, R, n_cells, n_ref, n_slices, H, W);
    // gen_laminate(domain.grid, 0.5, n_slices, n_ref);
    // gen_laminate(domain.grid, 1, n_ref, H, n_slices);
    // dealii::GridGenerator::hyper_cube(domain.grid, 0.0, 1.0);
    // domain.grid .refine_global (2);
    
    // arr<dbl, 3> size = {W, 1.0, H};
    cdbl P = 1.0;
    auto move = solve_cpeciment_extension_test (x, size, P, C, domain);
    // std::cout << move[x] << " " <<  move[y] << " " <<  move[z] << std::endl;
    // // std::cout << H/(move[z]*2.0) << std::endl;
    // std::cout << "Ez = " << H / (2.0 * move[z]) 
    //     << " ν_zx = " << -move[x] / move[z] << " ν_zy = " << -move[y] / move[z] << std::endl;
    // cdbl real_coef  = (R*R*M_PI*E2 + (1.0 - R*R*M_PI)*E1);
    // std::cout << real_coef << std::endl;
    // std::cout << std::abs(10.0/(move[z]*2.0) - real_coef) / real_coef * 100 << "%"  << std::endl;

    std::cout
        << "Ex = " << W / (2.0 * move[x]) 
        << " ν_xz = " << -(move[z] * W) / (move[x] * H)
        << " ν_xy = " << -(move[y] * W) / (move[x] * 1) << std::endl;

        // dealii::Triangulation<2> tria2d;
        // // dealii::GridGenerator ::hyper_cube (tria2d, 0.0, 1.0);
        // gen_circle_in_square_true_ordered(tria2d, R, n_ref);
        // // dealii::GridGenerator ::hyper_rectangle(tria2d, dealii::Point<2>(0.0, 0.0), dealii::Point<2>(W, 1.0));
        // // tria2d .refine_global (1);
        // std::cout << tria2d.n_vertices() << std::endl;
        // std::cout << tria2d.n_used_vertices() << std::endl;
        // auto v = tria2d.get_vertices();
        //
        // vec<dealii::CellData<2>> c (tria2d.n_active_cells(), dealii::CellData<2>());
        // {
        //     st cn = 0;
        //     auto cell = tria2d.begin_active();
        //     auto end_cell = tria2d.end();
        //     for (; cell != end_cell; ++cell)
        //     {
        //         for (st i = 0; i < 4; ++i)
        //         {
        //             c[cn].vertices[i] = cell->vertex_index(i);
        //         };
        //         ++cn;
        //     };
        // };
        // cst c_size = c.size();
        // for (st i = 0; i < c_size; ++i)
        // {
        //     dealii::CellData<2> new_c;
        //     for (st j = 0; j < 4; ++j)
        //     {
        //         auto new_v = v[c[i].vertices[j]];
        //         new_v(0) += 1.0;
        //         st indx = 0;
        //         bool new_v_exist = false;
        //         for (st k = 0; k < v.size(); ++k)
        //         {
        //             if (v[k].distance(new_v) < 1.0e-8)
        //             {
        //                 indx = k;
        //                 new_v_exist = true;
        //                 break;
        //             };
        //         };
        //         if (not new_v_exist)
        //         {
        //             v .push_back(new_v);
        //             indx = v.size()-1;
        //         };
        //         new_c.vertices[j] = indx;
        //     };
        //     new_c .material_id = c[i] .material_id;
        //     c .push_back(new_c);
        // };
        // std::cout << v.size() << std::endl;
        // std::cout << c.size() << std::endl;
        // dealii::Triangulation<2> tria2d_new;
        // tria2d_new.create_triangulation(v, c, dealii::SubCellData());
        // {
        //     std::ofstream out ("grid-speciment.eps");
        //     dealii::GridOut grid_out;
        //     grid_out.write_eps (tria2d_new, out);
        // };


    return 0;
}
