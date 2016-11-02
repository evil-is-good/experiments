#ifndef GRID_GENERATOR_gn584e48hf47838
#define GRID_GENERATOR_gn584e48hf47838 1
#include "feminist/calculation_core/src/blocks/general/point/point.h"

namespace CGALGridGenerator
{
    extern void set_grid(
            dealii::Triangulation< 2 >&,
            vec<prmt::Point<2>>,
            vec<vec<prmt::Point<2>>>,
            vec<st>);
    extern void set_grid_with_constraints(
            dealii::Triangulation< 2 >&,
            vec<prmt::Point<2>>,
            vec<vec<prmt::Point<2>>>,
            vec<st>);
    extern int QWERTY;
};

namespace GridGenerator
{
    void give_line_without_end_point(
            vec<prmt::Point<2>> &curve,
            cst num_points,
            prmt::Point<2> first,
            prmt::Point<2> second)
    {
        dbl dx = (second.x() - first.x()) / num_points;
        dbl dy = (second.y() - first.y()) / num_points;
        dbl x = first.x();
        dbl y = first.y();
        FOR_I(0, num_points - 0)
        {
            curve .push_back (prmt::Point<2>(x, y)); 
            x += dx;
            y += dy;
        };
    };

    void give_rectangle(
            vec<prmt::Point<2>> &curve,
            cst num_points_on_edge,
            prmt::Point<2> first,
            prmt::Point<2> second)
    {
        give_line_without_end_point(curve, num_points_on_edge,
                first,
                prmt::Point<2>(first.x(), second.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(first.x(), second.y()),
                second);

        give_line_without_end_point(curve, num_points_on_edge,
                second,
                prmt::Point<2>(second.x(), first.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(second.x(), first.y()),
                first);

    };

    void give_rectangle_with_border_condition(
            vec<prmt::Point<2>> &curve,
            vec<st> &type_edge,
            const arr<st, 4> type_border,
            cst num_points_on_edge,
            const prmt::Point<2> first,
            const prmt::Point<2> second)
    {
        give_line_without_end_point(curve, num_points_on_edge,
                first,
                prmt::Point<2>(first.x(), second.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(first.x(), second.y()),
                second);

        give_line_without_end_point(curve, num_points_on_edge,
                second,
                prmt::Point<2>(second.x(), first.y()));

        give_line_without_end_point(curve, num_points_on_edge,
                prmt::Point<2>(second.x(), first.y()),
                first);

        cst n_edge_on_border = curve.size() / 4;
        // printf("type %d\n", n_edge_on_border);
        type_edge.resize(curve.size());

        FOR(i, 0, 4)
            FOR(j, 0 + n_edge_on_border * i, n_edge_on_border + n_edge_on_border * i)
            type_edge[j] = type_border[i];
    };

    void give_circ(
            vec<prmt::Point<2>> &curve,
            cst num_points_on_tip,
            cdbl radius,
            prmt::Point<2> center)
    {
        cdbl angle_step_rad = 2.0 * M_PI / num_points_on_tip;
        for (
                dbl angle_rad = M_PI / 2.0; 
                std::abs(angle_rad - 5.0 * (M_PI / 2.0)) > 1.e-8; 
                angle_rad += angle_step_rad
            )
        {
            dbl X = radius * cos(angle_rad) + center.x();
            dbl Y = radius * sin(angle_rad) + center.y();
            curve .push_back (prmt::Point<2>(X, Y)); 
        };
    };

    void set_circles_in_rectangular_cgal(dealii::Triangulation<2> &triangulation,
            const arr<dbl, 3>& size,
            vec<dealii::Point<2>>& center,
            const double radius, cst n_points_on_includ_border)
    {
        vec<prmt::Point<2>> border;
        vec<st> type_border;
        give_rectangle_with_border_condition(
                border,
                type_border,
                arr<st, 4>{1,3,2,4},
                10,
                prmt::Point<2>(0.0, 0.0), prmt::Point<2>(size[0], size[1]));
        vec<vec<prmt::Point<2>>> inclusion(center.size());
        for (st i = 0; i < center.size(); ++i)
        {
            give_circ(inclusion[i], n_points_on_includ_border, radius, prmt::Point<2>(center[i]));
        };
        CGALGridGenerator::set_grid(triangulation, border, inclusion, type_border);
        {
            std::ofstream out ("grid-cgal.eps");
            dealii::GridOut grid_out;
            grid_out.write_eps (triangulation, out);
        };

        {
            auto cell = triangulation .begin_active();
            auto end_cell = triangulation .end();
            for (; cell != end_cell; ++cell)
            {
                cell->set_material_id(0);
                for (st i = 0; i < center.size(); ++i)
                {
                    if (center[i].distance(cell->center()) < radius)
                    {
                        cell->set_material_id(1);
                    };
                };
            };
        };
    };

    void set_cylinder_cgal(dealii::Triangulation<3> &triangulation, vec<dealii::Point<2>>& center, 
            const double radius, cst n_points_on_includ_border, cst n_slices)
    {
        dealii::Triangulation<2> tria2d;

        arr<dbl, 3> size = {1.0, 1.0, 1.0};
        set_circles_in_rectangular_cgal (tria2d, size, center, radius, n_points_on_includ_border); 

        GridGenerator::extrude_triangulation (tria2d, n_slices, size[2], triangulation);
    };

    void set_cylinder_in_rectangular_cgal(dealii::Triangulation<3> &triangulation,
            const arr<dbl, 3>& size,
            vec<dealii::Point<2>>& center, 
            const double radius, cst n_points_on_includ_border, cst n_slices)
    {
        dealii::Triangulation<2> tria2d;

        set_circles_in_rectangular_cgal (tria2d, size, center, radius, n_points_on_includ_border); 

        GridGenerator::extrude_triangulation (tria2d, n_slices, size[2], triangulation);
    };

    void mesh_tiling_2d (
            dealii::Triangulation<2>& tria,
            const dealii::Triangulation<2>& tile, 
            const arr<dbl, 2>& tile_size, cst n_tile_x, cst n_tile_y)
    {
        enum {x, y, z};
        std::cout << tile.n_vertices() << std::endl;

        cst n_tiles = n_tile_x * n_tile_y;
        cst n_tile_cells = tile.n_active_cells();

        auto v_tile = tile.get_vertices();
        vec<dealii::CellData<2>> c_tile (tile.n_active_cells(), dealii::CellData<2>());
        {
            auto cell = tile.begin_active();
            auto end_cell = tile.end();
            for (st cn = 0; cell != end_cell; ++cell)
            {
                for (st i = 0; i < 4; ++i)
                {
                    c_tile[cn].vertices[i] = cell->vertex_index(i);
                };
                ++cn;
            };
        };

        decltype(v_tile) v_tria;
        vec<dealii::CellData<2>> c_tria (tile.n_active_cells() * n_tiles, dealii::CellData<2>());

        {
            dbl shift_x = 0.0;
            st cn_tria = 0;
            for (st i = 0; i < n_tile_x ; ++i)
            {
                dbl shift_y = 0.0;
                for (st j = 0; j < n_tile_y; ++j)
                {
                    for (st cn = 0; cn < c_tile.size(); ++cn)
                    {
                        dealii::CellData<2> new_c;
                        for (st v = 0; v < 4; ++v)
                        {
                            auto new_v = v_tile[c_tile[cn].vertices[v]];

                            new_v(x) += shift_x;
                            new_v(y) += shift_y;

                            st indx = 0;
                            bool new_v_exist = false;
                            for (st k = 0; k < v_tria.size(); ++k)
                            {
                                if (v_tria[k].distance(new_v) < 1.0e-8)
                                {
                                    indx = k;
                                    new_v_exist = true;
                                    break;
                                };
                            };
                            if (not new_v_exist)
                            {
                                v_tria .push_back(new_v);
                                indx = v_tria.size()-1;
                            };
                            new_c.vertices[v] = indx;
                        };
                        new_c .material_id = c_tile[cn] .material_id;
                        c_tria[cn_tria] = new_c;
                        ++cn_tria;
                    };
                    shift_y += tile_size[y];
                };
                shift_x += tile_size[x];
            };
        };
        std::cout << "11111111111111" << std::endl;

        std::cout << v_tria.size() << std::endl;
        std::cout << c_tria.size() << std::endl;
        tria.create_triangulation(v_tria, c_tria, dealii::SubCellData());
        {
            std::ofstream out ("grid-speciment.eps");
            dealii::GridOut grid_out;
            grid_out.write_eps (tria, out);
        };
    };
};
#endif
