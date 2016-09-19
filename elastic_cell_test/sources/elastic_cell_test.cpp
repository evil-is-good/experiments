#include "pch.h"

#include "../../feminist/calculation_core/src/blocks/general/grid_generator/grid_generator.h"
#include "../../bloks/src/grid_generator/grid_generator.h"

ATools::FourthOrderTensor unphysical_to_physicaly (
        ATools::FourthOrderTensor &unphys)
{
    enum {x, y, z};
    ATools::FourthOrderTensor res;
    for (st i = 0; i < 3; ++i)
    {
        for (st j = 0; j < 3; ++j)
        {
            for (st k = 0; k < 3; ++k)
            {
                for (st l = 0; l < 3; ++l)
                {
                    res[i][j][k][l] = 0.0;
                };
            };
        };
    };

    double A = 
        unphys[x][x][x][x] * unphys[y][y][y][y] * unphys[z][z][z][z] - 
        unphys[y][y][z][z] * unphys[z][z][y][y] * unphys[x][x][x][x] +
        unphys[x][x][y][y] * unphys[y][y][z][z] * unphys[z][z][x][x] - 
        unphys[y][y][x][x] * unphys[x][x][y][y] * unphys[z][z][z][z] - 
        unphys[y][y][y][y] * unphys[x][x][z][z] * unphys[z][z][x][x] +
        unphys[y][y][x][x] * unphys[x][x][z][z] * unphys[z][z][y][y]; 

    for (uint8_t i = 0; i < 3; ++i)
    {
        int no_1 = (i + 1) % 3;
        int no_2 = (i + 2) % 3;

        for (uint8_t j = 0; j < 3; ++j)
        {
            int k = (j == no_1) ? no_2 : no_1;

            if (i == j)
                res[i][i][j][j] = A;
            else
                res[i][i][j][j] = 
                    (unphys[i][i][j][j] * unphys[k][k][k][k] -
                     unphys[i][i][k][k] * unphys[j][j][k][k]);

            res[i][i][j][j] /= 
                (unphys[no_1][no_1][no_1][no_1] * 
                 unphys[no_2][no_2][no_2][no_2] - 
                 unphys[no_1][no_1][no_2][no_2] * 
                 unphys[no_2][no_2][no_1][no_1]);
        };
    };

    return res;
};

void solve_approx_cell_elastic_problem (cdbl E_i, cdbl pua_i, cdbl E_m, cdbl pua_m,
        Domain<3> &domain)
{
        enum {x, y, z};

        dealii::FESystem<3,3> fe (dealii::FE_Q<3,3>(1), 3);
        domain.dof_init (fe);

        OnCell::SystemsLinearAlgebraicEquations<1> slae;
        OnCell::BlackOnWhiteSubstituter bows;

        LaplacianVector<3> element_matrix (domain.dof_handler.get_fe());
        element_matrix.C .resize (2);

        EPTools ::set_isotropic_elascity{yung : E_m, puasson : pua_m}(element_matrix.C[0]);
        EPTools ::set_isotropic_elascity{yung : E_i, puasson : pua_i}(element_matrix.C[1]);

        OnCell::prepare_system_equations_alternate<3, 3, 1> (slae, bows, domain);

        OnCell::Assembler::assemble_matrix<3> (slae.matrix, element_matrix, domain.dof_handler, bows);

        cst number_of_approx = 2; // Начиная с нулевой
        // arr<arr<i32, 3>, number_of_approx> approximations = {
        //     arr<i32, 3>({1, 0, 0}),
        //     arr<i32, 3>({0, 1, 0})};
        // arr<i32, 3>{2, 0, 0}};
        OnCell::ArrayWithAccessToVector<arr<arr<dbl, 3>, 3>> H(number_of_approx+1);
        OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> cell_func (number_of_approx);
        OnCell::ArrayWithAccessToVector<arr<dealii::Vector<dbl>, 3>> N_func (number_of_approx);
        OnCell::ArrayWithAccessToVector<arr<arr<dealii::Vector<dbl>, 3>, 3>> cell_stress (number_of_approx);
        OnCell::ArrayWithAccessToVector<arr<arr<dealii::Vector<dbl>, 3>, 3>> cell_deform (number_of_approx);
        OnCell::ArrayWithAccessToVector<arr<arr<arr<dbl, 3>, 3>, 3>> true_meta_coef (number_of_approx);
        for (auto &&a : H.content)
            for (auto &&b : a)
                for (auto &&c : b)
                    for (auto &&d : c)
                        for (auto &&e : d)
                            e = 0.0;
        for (auto &&a : cell_func.content)
            for (auto &&b : a)
                for (auto &&c : b)
                    for (auto &&d : c)
                        d .reinit (slae.solution[0].size());

        // Нулевое приближение ячейковой функции для которой nu==aplha равно 1.0
        for (st i = 0; i < slae.solution[0].size(); ++i)
        {
            if ((i % 3) == x) cell_func[arr<i32, 3>{0, 0, 0}][x][i] = 1.0;
            if ((i % 3) == y) cell_func[arr<i32, 3>{0, 0, 0}][y][i] = 1.0;
            if ((i % 3) == z) cell_func[arr<i32, 3>{0, 0, 0}][z][i] = 1.0;
        };
        // cell_func[arr<i32, 3>{0, 0, 0}][x][x] = 1.0; // Нулевое приближение ячейковой функции для которой nu==aplha равно 1.0
        // cell_func[arr<i32, 3>{0, 0, 0}][y][y] = 1.0; // Нулевое приближение ячейковой функции для которой nu==aplha равно 1.0
        // cell_func[arr<i32, 3>{0, 0, 0}][z][z] = 1.0; // Нулевое приближение ячейковой функции для которой nu==aplha равно 1.0
        for (auto &&a : N_func.content)
            for (auto &&b : a)
                for (auto &&c : b)
                    for (auto &&d : c)
                        d .reinit (slae.solution[0].size());
        for (auto &&a : cell_stress.content)
            for (auto &&b : a)
                for (auto &&c : b)
                    for (auto &&d : c)
                        for (auto &&e : d)
                            e .reinit (slae.solution[0].size());
        for (auto &&a : cell_deform.content)
            for (auto &&b : a)
                for (auto &&c : b)
                    for (auto &&d : c)
                        for (auto &&e : d)
                            e .reinit (slae.solution[0].size());

        OnCell::MetaCoefficientElasticCalculator mc_calculator (
                domain.dof_handler, element_matrix.C, domain.dof_handler.get_fe());

        // for (auto &&approximation : approximations)
        // {
        // auto approximation = approximations[0];
        for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
        {
            for (i32 i = 0; i < approx_number+1; ++i)
            {
                for (i32 j = 0; j < approx_number+1; ++j)
                {
                    for (i32 k = 0; k < approx_number+1; ++k)
                    {
                        if ((i+j+k) == approx_number)
                        {
                            arr<i32, 3> approximation = {i, j, k};
                            for (st nu = 0; nu < 3; ++nu)
                            {
                                slae.solution[0] = 0.0;
                                slae.rhsv[0] = 0.0;

                                OnCell::SourceVectorApprox<3> element_rhsv (approximation, nu,
                                        element_matrix.C, 
                                        H,
                                        cell_func,
                                        // &psi_func,
                                        domain.dof_handler.get_fe());
                                OnCell::Assembler::assemble_rhsv<3> (slae.rhsv[0], element_rhsv, domain.dof_handler, bows);

                                printf("problem %d %d %d %ld\n", i, j, k, nu);
                                printf("Integ %f\n", element_rhsv.tmp);
                                dealii::SolverControl solver_control (500000, 1e-12);
                                dealii::SolverCG<> solver (solver_control);
                                solver.solve (
                                        slae.matrix,
                                        slae.solution[0],
                                        slae.rhsv[0]
                                        ,dealii::PreconditionIdentity()
                                        );
                                FOR(i, 0, slae.solution[0].size())
                                    slae.solution[0][i] = slae.solution[0][bows.subst (i)];
                                FOR(i, 0, slae.rhsv[0].size())
                                    slae.rhsv[0][i] = slae.rhsv[0][bows.subst (i)];

                                cell_func[approximation][nu] = slae.solution[0];
                                // N_func[approximation][nu] = slae.rhsv[0];
                            };
                        };
                    };
                };
            };
            puts("!!!");
            for (i32 i = 0; i < approx_number+2; ++i)
            {
                for (i32 j = 0; j < approx_number+2; ++j)
                {
                    for (i32 k = 0; k < approx_number+2; ++k)
                    {
                        if ((i+j+k) == approx_number+1)
                        {
                            arr<i32, 3> approximation = {i, j, k};
                            for (st nu = 0; nu < 3; ++nu)
                            {
                                auto res = mc_calculator .calculate (
                                        approximation, nu,
                                        domain.dof_handler, cell_func);
                                H[approximation][nu][x] = res[x]; //E_x_a[0]_nu_a[1]
                                H[approximation][nu][y] = res[y]; 
                                H[approximation][nu][z] = res[z]; 
                                // printf("H k=(%d, %d, %d) nu=%ld %f %f %f\n", i, j, k, nu, 
                                //         H[approximation][nu][x],
                                //         H[approximation][nu][y],
                                //         H[approximation][nu][z]
                                //       );
                            };
                        };
                    };
                };
            };
        };
        puts("!!!!!");


        printf("\n");
        {
            arr<str, 3> ort = {"x", "y", "z"};
            arr<str, 3> aprx = {"0", "1", "2"};
            OnCell::StressCalculatorBD stress_calculator (
                    domain.dof_handler, element_matrix.C, domain.dof_handler.get_fe());
            OnCell::DeformCalculatorBD deform_calculator (
                    domain.dof_handler, domain.dof_handler.get_fe());
            for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
            {
                for (i32 i = 0; i < approx_number+1; ++i)
                {
                    for (i32 j = 0; j < approx_number+1; ++j)
                    {
                        for (i32 k = 0; k < approx_number+1; ++k)
                        {
                            if ((i+j+k) == approx_number)
                            {
                                arr<i32, 3> approximation = {i, j, k};
                                for (st nu = 0; nu < 3; ++nu)
                                {
                                    for (st alpha = 0; alpha < 3; ++alpha)
                                    {
                                        dealii::Vector<dbl> stress(domain.dof_handler.n_dofs());
                                        dealii::Vector<dbl> deform(domain.dof_handler.n_dofs());
                                        for (st beta = 0; beta < 3; ++beta)
                                        {
                                            // stress_calculator .calculate (
                                            //         approximation, nu, alpha, beta,
                                            //         domain.dof_handler, cell_func, stress);
                                            stress_calculator .calculate (
                                                    approximation, nu, alpha, beta,
                                                    cell_func, stress);
                                            deform_calculator .calculate (
                                                    approximation, nu, alpha, beta,
                                                    cell_func, deform);
                                            true_meta_coef[approximation][nu][alpha][beta] =
                                                OnCell::calculate_meta_coefficients_3d_elastic_from_stress (
                                                        domain.dof_handler, stress, beta);
                                            // printf("meta k=(%d, %d, %d) nu=%ld alpha=%ld beta=%ld %f\n", 
                                            //         i, j, k, nu, alpha, beta, true_meta_coef[approximation][nu][alpha][beta]);
                                        };
                                        cell_stress[approximation][nu][alpha] = stress;
                                        cell_deform[approximation][nu][alpha] = deform;
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };

        ATools::FourthOrderTensor meta_coef;
        for (st i = 0; i < 3; ++i)
        {
            for (st j = 0; j < 3; ++j)
            {
                for (st k = 0; k < 3; ++k)
                {
                    meta_coef[j][k][i][x] = true_meta_coef[arr<i32, 3>{1, 0, 0}][i][j][k];
                    meta_coef[j][k][i][y] = true_meta_coef[arr<i32, 3>{0, 1, 0}][i][j][k];
                    meta_coef[j][k][i][z] = true_meta_coef[arr<i32, 3>{0, 0, 1}][i][j][k];
                };
            };
        };
        //
        for (size_t i = 0; i < 9; ++i)
        {
            uint8_t im = i / (2 + 1);
            uint8_t in = i % (2 + 1);

            for (size_t j = 0; j < 9; ++j)
            {
                uint8_t jm = j / (2 + 1);
                uint8_t jn = j % (2 + 1);

                if (std::abs(meta_coef[im][in][jm][jn]) > 0.0000001)
                    printf("\x1B[31m%f\x1B[0m   ", 
                            meta_coef[im][in][jm][jn]);
                else
                    printf("%f   ", 
                            meta_coef[im][in][jm][jn]);
            };
            for (size_t i = 0; i < 2; ++i)
                printf("\n");
        };
        {
            auto newcoef = unphysical_to_physicaly (meta_coef);
            printf("%f %f %f %f %f %f %f %f %f %f %f\n", 
                    newcoef[0][0][0][0],
                    newcoef[0][0][1][1],
                    newcoef[0][0][2][2],
                    newcoef[1][1][0][0],
                    newcoef[1][1][1][1],
                    newcoef[1][1][2][2],
                    newcoef[2][2][0][0],
                    newcoef[2][2][1][1],
                    newcoef[2][2][2][2],
                    meta_coef[0][1][0][1],
                    meta_coef[0][2][0][2]
                  );
            // FILE* F;
            // F = fopen("coefs_ball_true.gpd", "a");
            // fprintf(F, "%f %f %f %f %f %f\n",
            //         // R*R*M_PI, 
            //         newcoef[0][0][0][0],
            //         newcoef[0][0][1][1],
            //         newcoef[0][0][2][2],
            //         newcoef[2][2][2][2],
            //         meta_coef[0][1][0][1],
            //         meta_coef[0][2][0][2]
            //         );
            // fclose(F);
        };

        {
            std::ofstream out ("cell/meta_coef.bin", std::ios::out | std::ios::binary);
            out.write ((char *) &meta_coef, sizeof meta_coef);
            out.close ();
        };

        arr<str, 3> ort = {"x", "y", "z"};
        arr<str, 3> aprx = {"0", "1", "2"};
        for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
        {
            for (i32 i = 0; i < approx_number+1; ++i)
            {
                for (i32 j = 0; j < approx_number+1; ++j)
                {
                    for (i32 k = 0; k < approx_number+1; ++k)
                    {
                        if ((i+j+k) == approx_number)
                        {
                            arr<i32, 3> approximation = {i, j, k};
                            for (st nu = 0; nu < 3; ++nu)
                            {
                                for (st alpha = 0; alpha < 3; ++alpha)
                                {
                                    str name = aprx[i]+str("_")+aprx[j]+str("_")+aprx[k]+str("_")+ort[nu]+str("_")+ort[alpha];
                                    {
                                        std::ofstream out ("cell/stress_"+name+".bin", std::ios::out | std::ios::binary);
                                        for (st i = 0; i < slae.solution[0].size(); ++i)
                                        {
                                            out.write ((char *) &(cell_stress[approximation][nu][alpha][i]), sizeof(dbl));
                                        };
                                        out.close ();
                                        EPTools ::print_move_slice (cell_stress[arr<i32,3>{i,j,k}][nu][alpha], domain.dof_handler, 
                                                "cell/stress_"+name+".gpd", z, 0.5);
                                    };
                                    {
                                        std::ofstream out ("cell/deform_"+name+".bin", std::ios::out | std::ios::binary);
                                        for (st i = 0; i < slae.solution[0].size(); ++i)
                                        {
                                            out.write ((char *) &(cell_deform[approximation][nu][alpha][i]), sizeof(dbl));
                                        };
                                        out.close ();
                                        EPTools ::print_move_slice (cell_deform[arr<i32,3>{i,j,k}][nu][alpha], domain.dof_handler, 
                                                "cell/deform_"+name+".gpd", z, 0.5);
                                    };
                                };
                            };
                        };
                    };
                };
            };
        };

        for (st approx_number = 1; approx_number < number_of_approx; ++approx_number)
        {
            for (i32 i = 0; i < approx_number+1; ++i)
            {
                for (i32 j = 0; j < approx_number+1; ++j)
                {
                    for (i32 k = 0; k < approx_number+1; ++k)
                    {
                        if ((i+j+k) == approx_number)
                        {
                            arr<i32, 3> approximation = {i, j, k};
                            for (st nu = 0; nu < 3; ++nu)
                            {
                                str name = ort[i]+str("_")+ort[j]+str("_")+ort[k]+str("_")+ort[nu];
                                {
                                    std::ofstream out ("cell/solution_"+name+".bin", std::ios::out | std::ios::binary);
                                    for (st i = 0; i < slae.solution[0].size(); ++i)
                                    {
                                        out.write ((char *) &(cell_func[approximation][nu][i]), sizeof(dbl));
                                    };
                                    out.close ();
                                };
                            };
                        };
                    };
                };
            };
        };

        {
            std::ofstream out ("cell/solution_on_cell_size.bin", std::ios::out | std::ios::binary);
            auto size = slae.solution[0].size();
            out.write ((char *) &size, sizeof size);
            out.close ();
        };



        OnCell::ArrayWithAccessToVector<arr<str, 3>> file_name (number_of_approx);
        file_name[arr<i32, 3>{1, 0, 0}][x] = "move_slice_approx_xx.gpd";
        file_name[arr<i32, 3>{1, 0, 0}][y] = "move_slice_approx_xy.gpd";
        file_name[arr<i32, 3>{1, 0, 0}][z] = "move_slice_approx_xz.gpd";
        file_name[arr<i32, 3>{0, 1, 0}][x] = "move_slice_approx_yx.gpd";
        file_name[arr<i32, 3>{0, 1, 0}][y] = "move_slice_approx_yy.gpd";
        file_name[arr<i32, 3>{0, 1, 0}][z] = "move_slice_approx_yz.gpd";
        EPTools ::print_coor_bin<3> (domain.dof_handler, "cell/coor_cell.bin");

        EPTools ::print_move_slice (cell_func[arr<i32, 3>{1, 0, 0}][x], domain.dof_handler, 
                file_name[arr<i32, 3>{1, 0, 0}][x], z, 0.5);
        EPTools ::print_move_slice (cell_func[arr<i32, 3>{1, 0, 0}][y], domain.dof_handler, 
                file_name[arr<i32, 3>{1, 0, 0}][y], z, 0.5);
        EPTools ::print_move_slice (cell_func[arr<i32, 3>{1, 0, 0}][z], domain.dof_handler, 
                file_name[arr<i32, 3>{1, 0, 0}][z], z, 0.5);
        EPTools ::print_move_slice (cell_func[arr<i32, 3>{0, 1, 0}][x], domain.dof_handler, 
                file_name[arr<i32, 3>{0, 1, 0}][x], z, 0.5);
        EPTools ::print_move_slice (cell_func[arr<i32, 3>{0, 1, 0}][y], domain.dof_handler, 
                file_name[arr<i32, 3>{0, 1, 0}][y], z, 0.5);
        EPTools ::print_move_slice (cell_func[arr<i32, 3>{0, 1, 0}][z], domain.dof_handler, 
                file_name[arr<i32, 3>{0, 1, 0}][z], z, 0.5);
        // EPTools ::print_move_slice (cell_func[arr<i32, 3>{2, 0, 0}][x], domain.dof_handler, 
        //         "move_slice_approx_2x_x.gpd", z, 0.5);
        // EPTools ::print_move_slice (cell_func[arr<i32, 3>{0, 2, 0}][x], domain.dof_handler, 
        //         "move_slice_approx_2y_x.gpd", z, 0.5);

        EPTools ::print_move_slice (cell_stress[arr<i32, 3>{1, 0, 0}][x][x], domain.dof_handler, 
                "stress_slice_xxx.gpd", z, 0.5);
        EPTools ::print_move_slice (cell_stress[arr<i32, 3>{0, 1, 0}][y][x], domain.dof_handler, 
                "stress_slice_yyx.gpd", z, 0.5);
};

void set_ball_true(dealii::Triangulation< 3 > &triangulation, 
        const double radius, const size_t n_refine)
{
    dealii::Point<3> center (0.5, 0.5, 0.5);
    dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
    triangulation .refine_global (n_refine);

    cdbl cell_size = 1.0 / std::pow(2.0, n_refine);
    cdbl max_dist = std::sqrt(std::pow(cell_size, 2.0) * 3.0) / 2.0;

    // auto face = tria2d.begin_active_face();
    // auto endf = tria2d.end_face();
    // for (; face != endf; ++face)
    {
        auto cell = triangulation .begin_active();
        auto end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            for (st i = 0; i < dealii::GeometryInfo<3>::vertices_per_cell; ++i)
            {
                auto p = cell->vertex(i);
                cdbl r = center.distance(p); 
                cdbl dist = std::abs(r - radius);
                // cdbl dist = std::abs(center.distance(p) - R); 
                if (dist < max_dist)
                {
                    cell->vertex(i) -= center;
                    cell->vertex(i) *= radius/r;
                    cell->vertex(i) += center;
                };
            };
        };
    };

    {
        auto cell = triangulation .begin_active();
        auto end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            if (center.distance(cell->center()) < radius)
            {
                cell->set_material_id(1);
            }
            else
                cell->set_material_id(0);
        };
    };
};

int main()
{
    enum {x, y, z};
    // CGALGridGenerator::QWERTY = 10;
    Domain<3> domain;
    cst n_ref = 4;
    // cdbl R = 0.415;//sqrt(0.7 / M_PI);//0.45;
    cdbl R = sqrt(0.7 / M_PI) / 4.0;//0.45;
    cst n_slices = 5;
    cst n_p = 16;
    // GridGenerator::gen_cylinder_in_cube_true_ordered(domain.grid, R, n_ref, n_slices);
    vec<dealii::Point<2>> center;
    center .push_back(dealii::Point<2>(0.5, 0.5));
    arr<dbl, 3> size = {1.0, 1.0, 1.0};
    GridGenerator::set_cylinder_in_rectangular_cgal(domain.grid, size, center, R, n_p, n_slices);
    // set_ball_true (domain.grid, R, n_ref);
    cdbl Em = 0.6;
    cdbl Ei = 60.0;
    cdbl pm = 0.35;
    cdbl pi = 0.2;
    // cdbl Em = 0.6;
    // cdbl Ei = 0.6;
    // cdbl pm = 0.25;
    // cdbl pi = 0.25;
    // cdbl Em = 0.3;
    // cdbl Ei = 60.0;
    // cdbl pm = 0.45;
    // cdbl pi = 0.2;
    solve_approx_cell_elastic_problem(Ei, pi, Em, pm, domain);

    return 0;
}
