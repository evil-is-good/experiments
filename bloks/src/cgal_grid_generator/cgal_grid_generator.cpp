/*
 * =====================================================================================
 *
 *       Filename:  grid.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01.07.2013 17:21:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef CGAL_GRID_GENERATOR_NF7E34WUDN
#define CGAL_GRID_GENERATOR_NF7E34WUDN 1

// #include "cgal/cgal.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
// #include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polygon_2.h>
// #include <CGAL/algorithm.h>

#include <stdlib.h>
#include <stdint.h>
#include <fstream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>

#include "../../../feminist/calculation_core/src/blocks/general/point/point.h"

namespace CGALGridGenerator
{
    struct FaceInfo2
    {
        FaceInfo2() : nesting_level(0) {};
        int nesting_level;

        bool in_domain(){ 
            return nesting_level%2 == 1;
        }
    };

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Triangulation_vertex_base_2<K> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
    // typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
    // typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
    // typedef CGAL::Triangulation_face_base_with_info_2 <FaceInfo2, K>    Fbb;
    // typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

    typedef CDT::Vertex_handle Vertex_handle;
    // typedef CDT::Point Point;

    void convert_to_dealii_format(CDT &cdt, dealii::Triangulation<2> &triangulation)
    {
        enum {d0, d1, d2}; // 0D, 1D, 2D

        std::vector<dealii::Point<2> > vertex_in_grid;
        std::vector<dealii::CellData<2> > cell_in_grid;

        for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        {
            CDT::Triangle trg = cdt.triangle(fit);

            // location of points on the triangle
            //
            //                    2
            //                   / \
            //                  /   \
            //                 /     \
            //                /       \
            //             5 /         \ 4
            //              /\         /\
            //             /  \       /  \
            //            /    \     /    \
            //           /      \   /      \
            //          /         6         \
            //         /          |          \
            //        /           |           \
            //       /            |            \
            //      0___________________________1
            //                    3
            //

            dealii::Point<2> p[7] =
            {
                dealii::Point<2>(trg[0].x(), trg[0].y()),
                dealii::Point<2>(trg[1].x(), trg[1].y()),
                dealii::Point<2>(trg[2].x(), trg[2].y()),

                (p[0] + p[1]) / 2.0,
                (p[1] + p[2]) / 2.0,
                (p[2] + p[0]) / 2.0,

                (p[0] + p[1] + p[2]) / 3.0
            };

            struct IndexPoint { u32 index; bool point_already_exist = false; };

            IndexPoint index_point[7];

            FOR(i, 0, 7)
                FOR(j, 0, vertex_in_grid.size())
                if (p[i].distance (vertex_in_grid[j]) < 1.e-10)
                {
                    index_point[i].index = j;
                    index_point[i].point_already_exist = true;
                };

            FOR(i, 0, 7)
                if (not index_point[i].point_already_exist)
                {
                    vertex_in_grid .push_back (p[i]);
                    index_point[i].index = vertex_in_grid.size() - 1;
                    index_point[i].point_already_exist = true;
                };

            u8 mat_id = fit->is_in_domain() ? 0 : 1;

            auto add_cell = [&cell_in_grid, &index_point, mat_id] (arr<st, 4> &&indx) 
            {
                dealii::CellData<d2> cell;
                cell.vertices[0] = index_point[indx[0]].index;
                cell.vertices[1] = index_point[indx[1]].index;
                cell.vertices[2] = index_point[indx[2]].index;
                cell.vertices[3] = index_point[indx[3]].index;
                cell.material_id = mat_id;
                cell_in_grid .push_back (cell);
                // cell_in_grid .push_back (dealii::CellData<d2>({
                //         index_point[indx[0]].index,
                //         index_point[indx[1]].index,
                //         index_point[indx[2]].index,
                //         index_point[indx[3]].index,
                //         {.material_id = mat_id}}));
            };

            // add_cell (arr<st, 4>{0, 3, 6, 5});
            // add_cell (arr<st, 4>{1, 4, 6, 3});
            // add_cell (arr<st, 4>{2, 5, 6, 4});
            add_cell (arr<st, 4>{3, 0, 6, 5});
            add_cell (arr<st, 4>{3, 1, 6, 4});
            add_cell (arr<st, 4>{6, 4, 5, 2});
        };
        {
            FILE *F;
            F = fopen("tria.gpd", "w");
            for (st i = 0; i < vertex_in_grid.size(); ++i)
            {
                fprintf(F, "%f %f 0.0\n", vertex_in_grid[i](0), vertex_in_grid[i](1));
            };
            fclose(F);
        };
        {
            FILE *F;
            F = fopen("9.gpd", "w");
            fprintf (F, "%f %f 0.0\n",
                    vertex_in_grid[cell_in_grid[9].vertices[0]](0),
                    vertex_in_grid[cell_in_grid[9].vertices[0]](1));
            fprintf (F, "%f %f 0.0\n",
                    vertex_in_grid[cell_in_grid[9].vertices[1]](0),
                    vertex_in_grid[cell_in_grid[9].vertices[1]](1));
            fprintf (F, "%f %f 0.0\n",
                    vertex_in_grid[cell_in_grid[9].vertices[2]](0),
                    vertex_in_grid[cell_in_grid[9].vertices[2]](1));
            fprintf (F, "%f %f 0.0\n",
                    vertex_in_grid[cell_in_grid[9].vertices[3]](0),
                    vertex_in_grid[cell_in_grid[9].vertices[3]](1));
            for (st i = 0; i < cell_in_grid.size(); ++i)
            {
                printf("%d %d %d %d %d\n", 
                        i,
                        cell_in_grid[i].vertices[0],
                        cell_in_grid[i].vertices[1],
                        cell_in_grid[i].vertices[2],
                        cell_in_grid[i].vertices[3]
                      );
            };
            fclose(F);
        };

        puts("111");
        // dealii::GridReordering<2> ::reorder_cells (cell_in_grid);
        puts("111");
        // triangulation .create_triangulation_compatibility (
        triangulation .create_triangulation (
                vertex_in_grid, cell_in_grid, dealii::SubCellData());
        puts("111");
    };

    // void make_edge_in_grid (
    //         const CDT &cdt, 
    //         const std::vector<dealii::Point<2>> &vertex_in_grid,
    //         const vec<prmt::Point<2>> &outer_border,
    //         const vec<st> &type_outer_border,
    //         dealii::SubCellData &edge_in_grid, 
    //         st &last_index)
    // {
    //     // FOR(i, 0, outer_border.size() - 1)
    //     // {
    //     //     vec<pair<st, CDT::Point>> point_on_edge;
    //     //     vec<st> index_on_edge;
    //     //     dbl l = outer_border[0].distance (outer_border[i+1]);
    //     //     FOR
    //     // };
    // };

    bool point_on_segment (
            const dealii::Point<2> &point, 
            const arr<dealii::Point<2>, 2> &segment)
    {
        cdbl x = point(0);
        cdbl y = point(1);
        cdbl x1 = segment[0](0);
        cdbl y1 = segment[0](1);
        cdbl x2 = segment[1](0);
        cdbl y2 = segment[1](1);

        cdbl l = segment[0].distance (segment[1]);

        if (abs(x * (y1 - y2) - y * (x1 - x2) + x1 * y2 - x2 * y1) < 1e-10)
            if (
                    ((l - segment[0].distance (point)) > -1e-10) and
                    ((l - segment[1].distance (point)) > -1e-10)
               )
                return true;

        return false;
    }

    void convert_to_dealii_format(
            const CDT &cdt, 
            const vec<CDT::Point> &border_on_cdt,
            dealii::Triangulation<2> &triangulation,
            const vec<prmt::Point<2>> &outer_border,
            const vec<st> &type_outer_border)
    {
        enum {d0, d1, d2}; // 0D, 1D, 2D

        std::vector<dealii::Point<2> > vertex_in_grid;
        std::vector<dealii::CellData<2> > cell_in_grid;
        dealii::SubCellData edge_in_grid;

        vec<st> index_boundary_point(outer_border.size());

        for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        {
            // fit->is_infinite();
            // fit->is_constrained(0);
            // if (fit->is_in_domain())
            //         {
            CDT::Triangle trg = cdt.triangle(fit);
            cdt.is_infinite (fit);

            // location of points on the triangle
            //
            //                    2
            //                   / \
            //                  /   \
            //                 /     \
            //                /       \
            //             5 /         \ 4
            //              /\         /\
            //             /  \       /  \
            //            /    \     /    \
            //           /      \   /      \
            //          /         6         \
            //         /          |          \
            //        /           |           \
            //       /            |            \
            //      0___________________________1
            //                    3
            //

            dealii::Point<2> p[7] =
            {
                dealii::Point<2>(trg[0].x(), trg[0].y()),
                dealii::Point<2>(trg[1].x(), trg[1].y()),
                dealii::Point<2>(trg[2].x(), trg[2].y()),

                (p[0] + p[1]) / 2.0,
                (p[1] + p[2]) / 2.0,
                (p[2] + p[0]) / 2.0,

                (p[0] + p[1] + p[2]) / 3.0
            };

            struct IndexPoint { u32 index; bool point_already_exist = false; };

            IndexPoint index_point[7];

            FOR(i, 0, 7)
                FOR(j, 0, vertex_in_grid.size())
                if (p[i].distance (vertex_in_grid[j]) < 1.e-10)
                {
                    index_point[i].index = j;
                    index_point[i].point_already_exist = true;
                };

            FOR(i, 0, 7)
                if (not index_point[i].point_already_exist)
                {
                    vertex_in_grid .push_back (p[i]);
                    index_point[i].index = vertex_in_grid.size() - 1;
                    index_point[i].point_already_exist = true;

                    FOR(j, 0, outer_border.size())
                        if (outer_border[j] == p[i])
                            index_boundary_point[j] = index_point[i].index;
                };

            u8 mat_id = 0;//fit->is_in_domain() ? 0 : 1;

            auto add_cell = [&cell_in_grid, &index_point, mat_id] (arr<st, 4> &&indx) 
            {
                dealii::CellData<d2> cell;
                cell.vertices[0] = index_point[indx[0]].index;
                cell.vertices[1] = index_point[indx[1]].index;
                cell.vertices[2] = index_point[indx[2]].index;
                cell.vertices[3] = index_point[indx[3]].index;
                cell.material_id = mat_id;
                cell_in_grid .push_back (cell);
                // cell_in_grid .push_back (dealii::CellData<d2>{
                //         index_point[indx[0]].index,
                //         index_point[indx[1]].index,
                //         index_point[indx[2]].index,
                //         index_point[indx[3]].index,
                //         {.material_id = mat_id}});
            };

            // if (fit->is_in_domain())
            //     puts("in");
            // else
            //     puts("not in");
            add_cell (arr<st, 4>{0, 3, 6, 5});
            add_cell (arr<st, 4>{1, 4, 6, 3});
            add_cell (arr<st, 4>{2, 5, 6, 4});
            // };
        };

        // FOR (i, 0, index_boundary_point.size() - 1)
        //     edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
        //             index_boundary_point[i],
        //             index_boundary_point[i + 1],
        //             {.boundary_id = type_outer_border[i]}});

        // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
        //         index_boundary_point.back(),
        //         index_boundary_point.front(),
        //         {.boundary_id = type_outer_border.back()}});

        // create_file("test_grid.gpd");
        // append_in_file("test_grid.gpd", outer_border.size
        // FOR(i, 0, edeg_in_grid.boundary_lines.size())
        //     append_in_file("test_grid.gpd", edeg_in_grid.boundary_lines
        // st first_index = 0;
        // st last_index  = 0;
        // FOR(i, 0, vertex_in_grid.size())
        //     if (vertex_in_grid[i].distance(
        //                 dealii::Point<2>(outer_border[0].x(), outer_border[0].y()))
        //             < 1e-10)
        //         first_index = i;
        // printf("first_id %d\n", first_index);

        // make_edge_in_grid (
        //         cdt, 
        //         vertex_in_grid, 
        //         outer_border, 
        //         type_outer_border, 
        //         edge_in_grid, 
        //         last_index);

        FOR(i, 0, border_on_cdt.size() - 1)
        {
            st index_1 = 0xFFFFFFFF;
            st index_2 = 0xFFFFFFFF;
            st index_3 = 0xFFFFFFFF;

            dealii::Point<2> first_point (
                    border_on_cdt[i].x(), 
                    border_on_cdt[i].y());

            dealii::Point<2> midl_point (
                    (border_on_cdt[i].x() + border_on_cdt[i + 1].x()) / 2.0,
                    (border_on_cdt[i].y() + border_on_cdt[i + 1].y()) / 2.0);

            dealii::Point<2> second_point (
                    border_on_cdt[i + 1].x(), 
                    border_on_cdt[i + 1].y());

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(first_point) < 1e-10)
                    index_1 = j;

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(midl_point) < 1e-10)
                    index_2 = j;

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(second_point) < 1e-10)
                    index_3 = j;

            st type_border = 0;
            bool fl = false;
            // printf("n %ld\n", outer_border.size());
            FOR(j, 0, outer_border.size())
            {
                const arr<dealii::Point<2>, 2> segment = {
                    outer_border[j], outer_border[(j + 1) % 4]};
                // printf("(%f,%f) (%f,%f) (%f,%f) (%f,%f)\n",
                //         segment[0](0), segment[0](1),
                //         segment[1](0), segment[1](1),
                //         vertex_in_grid[index_1](0),
                //         vertex_in_grid[index_1](1),
                //         vertex_in_grid[index_3](0),
                //         vertex_in_grid[index_3](1));

                if (
                        point_on_segment(vertex_in_grid[index_1], segment) and
                        point_on_segment(vertex_in_grid[index_3], segment))
                {
                    type_border = type_outer_border[j];
                    fl = true;
                    break;
                };
            };
            // printf("olol %ld %d\n", type_border, fl);


            {
                dealii::CellData<d1> cell;
                cell.vertices[0] = index_1;
                cell.vertices[1] = index_2;
                cell.boundary_id = type_border;
                edge_in_grid.boundary_lines .push_back (cell);
                // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
                //         index_1,
                //         index_2,
                //         {.boundary_id = type_border}});
            };

            {
                dealii::CellData<d1> cell;
                cell.vertices[0] = index_2;
                cell.vertices[1] = index_3;
                cell.boundary_id = type_border;
                edge_in_grid.boundary_lines .push_back (cell);
                // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
                //         index_2,
                //         index_3,
                //         {.boundary_id = type_border}});
            };

            if (
                    (index_1 == 0xFFFFFFFF) or
                    (index_2 == 0xFFFFFFFF) or
                    (index_3 == 0xFFFFFFFF) 
               )
                puts("\x1B[31mERRRRROOOOOOORRRR!!!!!!!111\x1B[0m");
        };

        dealii::GridReordering<2> ::reorder_cells (cell_in_grid);
        triangulation .create_triangulation_compatibility (
                vertex_in_grid, cell_in_grid, 
                // dealii::SubCellData());
                      edge_in_grid);
    };

    void convert_to_dealii_format_without_domains(
            const CDT &cdt, 
            const vec<CDT::Point> &border_on_cdt,
            dealii::Triangulation<2> &triangulation,
            const vec<prmt::Point<2>> &outer_border,
            const vec<st> &type_outer_border)
    {
        enum {d0, d1, d2}; // 0D, 1D, 2D

        std::vector<dealii::Point<2> > vertex_in_grid;
        std::vector<dealii::CellData<2> > cell_in_grid;
        dealii::SubCellData edge_in_grid;

        vec<st> index_boundary_point(outer_border.size());

        for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        {
            // fit->is_infinite();
            // fit->is_constrained(0);
            // if (fit->is_in_domain())
            //         {
            CDT::Triangle trg = cdt.triangle(fit);
            cdt.is_infinite (fit);

            // location of points on the triangle
            //
            //                    2
            //                   / \
            //                  /   \
            //                 /     \
            //                /       \
            //             5 /         \ 4
            //              /\         /\
            //             /  \       /  \
            //            /    \     /    \
            //           /      \   /      \
            //          /         6         \
            //         /          |          \
            //        /           |           \
            //       /            |            \
            //      0___________________________1
            //                    3
            //

            dealii::Point<2> p[7] =
            {
                dealii::Point<2>(trg[0].x(), trg[0].y()),
                dealii::Point<2>(trg[1].x(), trg[1].y()),
                dealii::Point<2>(trg[2].x(), trg[2].y()),

                (p[0] + p[1]) / 2.0,
                (p[1] + p[2]) / 2.0,
                (p[2] + p[0]) / 2.0,

                (p[0] + p[1] + p[2]) / 3.0
            };

            struct IndexPoint { u32 index; bool point_already_exist = false; };

            IndexPoint index_point[7];

            FOR(i, 0, 7)
                FOR(j, 0, vertex_in_grid.size())
                if (p[i].distance (vertex_in_grid[j]) < 1.e-10)
                {
                    index_point[i].index = j;
                    index_point[i].point_already_exist = true;
                };

            FOR(i, 0, 7)
                if (not index_point[i].point_already_exist)
                {
                    vertex_in_grid .push_back (p[i]);
                    index_point[i].index = vertex_in_grid.size() - 1;
                    index_point[i].point_already_exist = true;

                    FOR(j, 0, outer_border.size())
                        if (outer_border[j] == p[i])
                            index_boundary_point[j] = index_point[i].index;
                };

            u8 mat_id = 0;//fit->is_in_domain() ? 0 : 1;

            auto add_cell = [&cell_in_grid, &index_point, mat_id] (arr<st, 4> &&indx) 
            {
                dealii::CellData<d2> cell;
                cell.vertices[0] = index_point[indx[0]].index;
                cell.vertices[1] = index_point[indx[1]].index;
                cell.vertices[2] = index_point[indx[2]].index;
                cell.vertices[3] = index_point[indx[3]].index;
                cell.material_id = mat_id;
                cell_in_grid .push_back (cell);
                // cell_in_grid .push_back (dealii::CellData<d2>{
                //         index_point[indx[0]].index,
                //         index_point[indx[1]].index,
                //         index_point[indx[2]].index,
                //         index_point[indx[3]].index,
                //         {.material_id = mat_id}});
            };

            // if (fit->is_in_domain())
            //     puts("in");
            // else
            //     puts("not in");
            add_cell (arr<st, 4>{0, 3, 6, 5});
            add_cell (arr<st, 4>{1, 4, 6, 3});
            add_cell (arr<st, 4>{2, 5, 6, 4});
            // };
        };

        // FOR (i, 0, index_boundary_point.size() - 1)
        //     edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
        //             index_boundary_point[i],
        //             index_boundary_point[i + 1],
        //             {.boundary_id = type_outer_border[i]}});

        // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
        //         index_boundary_point.back(),
        //         index_boundary_point.front(),
        //         {.boundary_id = type_outer_border.back()}});

        // create_file("test_grid.gpd");
        // append_in_file("test_grid.gpd", outer_border.size
        // FOR(i, 0, edeg_in_grid.boundary_lines.size())
        //     append_in_file("test_grid.gpd", edeg_in_grid.boundary_lines
        // st first_index = 0;
        // st last_index  = 0;
        // FOR(i, 0, vertex_in_grid.size())
        //     if (vertex_in_grid[i].distance(
        //                 dealii::Point<2>(outer_border[0].x(), outer_border[0].y()))
        //             < 1e-10)
        //         first_index = i;
        // printf("first_id %d\n", first_index);

        // make_edge_in_grid (
        //         cdt, 
        //         vertex_in_grid, 
        //         outer_border, 
        //         type_outer_border, 
        //         edge_in_grid, 
        //         last_index);

        FOR(i, 0, border_on_cdt.size() - 1)
        {
            st index_1 = 0xFFFFFFFF;
            st index_2 = 0xFFFFFFFF;
            st index_3 = 0xFFFFFFFF;

            dealii::Point<2> first_point (
                    border_on_cdt[i].x(), 
                    border_on_cdt[i].y());

            dealii::Point<2> midl_point (
                    (border_on_cdt[i].x() + border_on_cdt[i + 1].x()) / 2.0,
                    (border_on_cdt[i].y() + border_on_cdt[i + 1].y()) / 2.0);

            dealii::Point<2> second_point (
                    border_on_cdt[i + 1].x(), 
                    border_on_cdt[i + 1].y());

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(first_point) < 1e-10)
                    index_1 = j;

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(midl_point) < 1e-10)
                    index_2 = j;

            FOR(j, 0, vertex_in_grid.size())
                if (vertex_in_grid[j].distance(second_point) < 1e-10)
                    index_3 = j;

            st type_border = 0;
            bool fl = false;
            printf("n %ld\n", outer_border.size());
            FOR(j, 0, outer_border.size())
            {
                const arr<dealii::Point<2>, 2> segment = {
                    outer_border[j], outer_border[(j + 1) % 4]};
                printf("(%f,%f) (%f,%f) (%f,%f) (%f,%f)\n",
                        segment[0](0), segment[0](1),
                        segment[1](0), segment[1](1),
                        vertex_in_grid[index_1](0),
                        vertex_in_grid[index_1](1),
                        vertex_in_grid[index_3](0),
                        vertex_in_grid[index_3](1));

                if (
                        point_on_segment(vertex_in_grid[index_1], segment) and
                        point_on_segment(vertex_in_grid[index_3], segment))
                {
                    type_border = type_outer_border[j];
                    fl = true;
                    break;
                };
            };
            printf("olol %ld %d\n", type_border, fl);

            {
                dealii::CellData<d1> cell;
                cell.vertices[0] = index_1;
                cell.vertices[1] = index_2;
                cell.boundary_id = type_border;
                edge_in_grid.boundary_lines .push_back (cell);
                // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
                //         index_1,
                //         index_2,
                //         {.boundary_id = type_border}});
            };

            {
                dealii::CellData<d1> cell;
                cell.vertices[0] = index_2;
                cell.vertices[1] = index_3;
                cell.boundary_id = type_border;
                edge_in_grid.boundary_lines .push_back (cell);
                // edge_in_grid.boundary_lines .push_back (dealii::CellData<d1>{
                //         index_2,
                //         index_3,
                //         {.boundary_id = type_border}});
            };

            if (
                    (index_1 == 0xFFFFFFFF) or
                    (index_2 == 0xFFFFFFFF) or
                    (index_3 == 0xFFFFFFFF) 
               )
                puts("\x1B[31mERRRRROOOOOOORRRR!!!!!!!111\x1B[0m");
        };

        dealii::GridReordering<2> ::reorder_cells (cell_in_grid);
        triangulation .create_triangulation_compatibility (
                vertex_in_grid, cell_in_grid, 
                // dealii::SubCellData());
                      edge_in_grid);
    };

    enum class cguc_type
    {
        Default
    };
    template <cguc_type type/*  = cguc_type::Default */>
        void create_grid_using_cgal(CDT &cdt,
                const vec<prmt::Point<2>> &outer_border,
                const vec<prmt::Point<2>> &inner_border)
        {};

    // void insert_polygon(CDT& cdt, const CGAL::Polygon_2<K>& polygon){
    //     if ( polygon.is_empty() ) return;
    //     CDT::Vertex_handle v_prev=cdt.insert(*CGAL::cpp0x::prev(polygon.vertices_end()));
    //     for (auto vit=polygon.vertices_begin();
    //             vit!=polygon.vertices_end();++vit)
    //     {
    //         CDT::Vertex_handle vh=cdt.insert(*vit);
    //         cdt.insert_constraint(vh,v_prev);
    //         v_prev=vh;
    //     }
    // }

    template <>
        void create_grid_using_cgal<cguc_type::Default> (CDT &cdt,
                const vec<prmt::Point<2>> &outer_border,
                const vec<prmt::Point<2>> &inner_border)
        {
            std::vector<std::vector<Vertex_handle> > vec_of_domains;

            auto add_border_to_domain = [&vec_of_domains, &cdt] (vec<prmt::Point<2>> border)
            {
                std::vector<Vertex_handle> vec_of_vertices;
                for(auto p : border)
                    vec_of_vertices.push_back(cdt.insert(CDT::Point(p.x(), p.y())));

                vec_of_domains.push_back(vec_of_vertices);
            };

            add_border_to_domain (inner_border);
            add_border_to_domain (outer_border);
            // 
            for(auto it = vec_of_domains.begin(); it != vec_of_domains.end(); ++it) 
            {
                // auto it = vec_of_domains.begin();
                // ++it;
                for(auto vit = it->begin() + 1; vit != it->end(); ++vit)
                {
                    cdt.insert_constraint( *(vit - 1), *vit );
                }
                cdt.insert_constraint( *( it->end() - 1), *( it->begin() ) );
            }

            std::list<CDT::Point> list_of_seeds;
            list_of_seeds.push_back(CDT::Point(0.5, 0.5));
            // // std::list<Point> list_of_seeds_1;
            // // list_of_seeds_1.push_back(Point(0.5, 0.5));
            // // list_of_seeds_1.push_back(Point(0.2, 0.2));
            // // std::list<Point> list_of_seeds_2;
            // // list_of_seeds_2.push_back(Point(0.5, 0.5));

            CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                    Criteria( 0.125, 0.1 ), true);
            CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                    Criteria(0.125, 0.1));

        };

    void extract_border (const CDT &cdt, vec<CDT::Point> &border_on_cdt)
    {
        // auto fc = cdt.incident_edges (cdt.infinite_vertex());
        CDT::Face_circulator fc = cdt.incident_faces (cdt.infinite_vertex());
        CDT::Face_circulator end_fc = fc;
        printf("2.1.2.1\n");
        FOR(i, 0, 9000)
        {
            st e = 3;
            printf("2.1.2.2\n");
            // if (fc->is_constrained(0)) e = 0;
            printf("2.1.2.3\n");
            // if (fc->is_constrained(1)) e = 1;
            // if (fc->is_constrained(2)) e = 2;
            for (st i = 0; i < 3; ++i)
            {
                if (fc->is_constrained(i))
                    e = i;
            };

            if (e < 3)
            {
                // append_in_file("border_grid_cgal.gpd", 
                //         std::to_string(fc->vertex((e+1) % 3)->point().x()) +
                //         " " +
                //         std::to_string(fc->vertex((e+1) % 3)->point().y()) +
                //         " 0.0\n"); 
                // append_in_file("border_grid_cgal.gpd", 
                //         std::to_string(fc->vertex((e+2) % 3)->point().x()) +
                //         " " +
                //         std::to_string(fc->vertex((e+2) % 3)->point().y()) +
                //         " 0.0\n"); 
                border_on_cdt .push_back (fc->vertex((e + 1) % 3)->point());
                // border_on_cdt .push_back (fc->vertex((e + 2) % 3)->point());
            };
            ++fc;
            if (fc == end_fc)
            {
                break;
            };
        };
        printf("2.1.2.4\n");
        printf("%d\n", border_on_cdt.size());
        border_on_cdt .push_back (border_on_cdt.front()); 
        printf("2.1.2.5\n");
        // for (auto i : border_on_cdt)
        //     append_in_file("border_grid_cgal.gpd", 
        //             std::to_string(i.x()) +
        //             " " +
        //             std::to_string(i.y()) +
        //             " 0.0\n"); 
    };

    void create_grid_using_cgal_1 (CDT &cdt,
            vec<CDT::Point> &border_on_cdt,
            const vec<prmt::Point<2>> &outer_border,
            const vec<prmt::Point<2>> &inner_border)
    {
        // CGAL::Polygon_2<K> polygon[2];

        // for (auto p : outer_border)
        //     polygon[0] .push_back (CDT::Point(p.x(), p.y()));
        // for (auto p : inner_border)
        //     polygon[1] .push_back (CDT::Point(p.x(), p.y()));

        // insert_polygon (cdt, polygon[0]);
        // insert_polygon (cdt, polygon[1]);
        std::vector<std::vector<Vertex_handle> > vec_of_domains;

        auto add_border_to_domain = [&vec_of_domains, &cdt] (vec<prmt::Point<2>> border)
        {
            std::vector<Vertex_handle> vec_of_vertices;
            for(auto p : border)
                vec_of_vertices.push_back(cdt.insert(CDT::Point(p.x(), p.y())));

            vec_of_domains.push_back(vec_of_vertices);
        };

        add_border_to_domain (inner_border);
        add_border_to_domain (outer_border);
        // 
        for(auto it = vec_of_domains.begin(); it != vec_of_domains.end(); ++it) 
        {
            // auto it = vec_of_domains.begin();
            // ++it;
            for(auto vit = it->begin() + 1; vit != it->end(); ++vit)
            {
                cdt.insert_constraint( *(vit - 1), *vit );
            }
            cdt.insert_constraint( *( it->end() - 1), *( it->begin() ) );
        }

        // for (auto e = cdt.finite_edges_begin(); e != cdt.finite_edges_end(); ++e)
        // {
        //     cdt.is_constrained(CDT::Edge(e));
        // };
        // 
        // create_file("border_grid_cgal.gpd");
        // for (auto i : border)
        // i.first->is_in_domain();
        // i.first->vertex(0);
        // append_in_file("border_grid_cgal.gpd", 
        //         std::to_string(i.first->vertex(2)->point().x()) +
        //         " " +
        //         std::to_string(i.first->vertex(2)->point().y()) +
        //         "\n"); 

        std::list<CDT::Point> list_of_seeds;
        list_of_seeds.push_back(CDT::Point(0.5, 0.5));
        // // std::list<Point> list_of_seeds_1;
        // // list_of_seeds_1.push_back(Point(0.5, 0.5));
        // // list_of_seeds_1.push_back(Point(0.2, 0.2));
        // // std::list<Point> list_of_seeds_2;
        // // list_of_seeds_2.push_back(Point(0.5, 0.5));

        CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                Criteria(/*  0.125, 0.1 */), true);
        CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                Criteria(/* 0.125, 0.1 */));
        // vec<CDT::Edge> border;
        // // CDT::Face_handle fh = cdt.infinite_face();
        // CDT::Vertex_handle v;
        // for(auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        // {
        //     st e = 3;
        //     if (fit->is_constrained(0)) e = 0;
        //     if (fit->is_constrained(1)) e = 1;
        //     if (fit->is_constrained(2)) e = 2;
        //     if (e < 3)
        //     {
        //         append_in_file("border_grid_cgal.gpd", 
        //                 std::to_string(fit->vertex((e+1) % 3)->point().x()) +
        //                 " " +
        //                 std::to_string(fit->vertex((e+1) % 3)->point().y()) +
        //                 " 0.0\n"); 
        //         append_in_file("border_grid_cgal.gpd", 
        //                 std::to_string(fit->vertex((e+2) % 3)->point().x()) +
        //                 " " +
        //                 std::to_string(fit->vertex((e+2) % 3)->point().y()) +
        //                 " 0.0\n"); 
        //         v = fit->vertex((e + 2) % 3);
        //         break;
        //     };
        // //     if (not fit->is_in_domain())
        // //     for(int i = 0; i < 3; i++)
        // // append_in_file("border_grid_cgal.gpd", 
        // //         std::to_string(fit->vertex(i)->point().x()) +
        // //         " " +
        // //         std::to_string(fit->vertex(i)->point().y()) +
        // //         " 0.0\n"); 
        //     // for(int i = 0; i < 3; i++)
        //     // {
        //     // CDT::Edge e(fh,i);
        //     // CDT::Face_handle n = fh->neighbor(i);
        //     // if (cdt.is_constrained(e)) border.push_back(e);
        //     // else fh = n;
        //     // };
        // };
        extract_border (cdt, border_on_cdt);
    };

    void create_grid_using_cgal_1 (CDT &cdt,
            vec<CDT::Point> &border_on_cdt,
            const vec<prmt::Point<2>> &outer_border,
            const vec<vec<prmt::Point<2>>> &inner_border)
    {
        std::vector<std::vector<Vertex_handle> > vec_of_domains;

        auto add_border_to_domain = [&vec_of_domains, &cdt] (vec<prmt::Point<2>> border)
        {
            std::vector<Vertex_handle> vec_of_vertices;
            for(auto p : border)
                vec_of_vertices.push_back(cdt.insert(CDT::Point(p.x(), p.y())));

            vec_of_domains.push_back(vec_of_vertices);
        };

        for (st i = 0; i < inner_border.size(); ++i)
        {
            add_border_to_domain (inner_border[i]);
        };
        add_border_to_domain (outer_border);
        // 
        for(auto it = vec_of_domains.begin(); it != vec_of_domains.end(); ++it) 
        {
            for(auto vit = it->begin() + 1; vit != it->end(); ++vit)
            {
                cdt.insert_constraint( *(vit - 1), *vit );
            }
            cdt.insert_constraint( *( it->end() - 1), *( it->begin() ) );
        }

        std::list<CDT::Point> list_of_seeds;
        list_of_seeds.push_back(CDT::Point(0.5, 0.5));

        CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                Criteria(/*  0.125, 0.1 */), true);
        CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                Criteria(/* 0.125, 0.1 */));
        extract_border (cdt, border_on_cdt);
    };

    void create_grid_using_cgal_2 (CDT &cdt,
            vec<CDT::Point> &border_on_cdt,
            const vec<prmt::Point<2>> &outer_border)
    {
        vec<vec<Vertex_handle>> vec_of_domains;

        auto add_border_to_domain = [&vec_of_domains, &cdt] (vec<prmt::Point<2>> border)
        {
            std::vector<Vertex_handle> vec_of_vertices;
            for(auto p : border)
                vec_of_vertices.push_back(cdt.insert(CDT::Point(p.x(), p.y())));

            vec_of_domains.push_back(vec_of_vertices);
        };

        auto add_constreins = [&vec_of_domains, &cdt] ()
        {
            for (auto vv : vec_of_domains)
            {
                for (st i = 0; i < vv.size() - 1; ++i)
                {
                    cdt.insert_constraint(vv[i], vv[i + 1]);
                };
                cdt.insert_constraint(vv.back(), vv.front());
            };
        };
        printf("2.1.1\n");

        add_border_to_domain (outer_border);
        add_constreins ();

        // CGAL::refine_Delaunay_mesh_2(cdt, Criteria( 0.125, 0.1 ));
        CGAL::refine_Delaunay_mesh_2(cdt, Criteria(/* 0.125, 0.1 */));
        printf("2.1.2\n");

        extract_border (cdt, border_on_cdt);
        printf("2.1.3\n");
    };

    void set_grid(
            dealii::Triangulation< 2 > &triangulation,
            vec<prmt::Point<2>> outer_border,
            vec<prmt::Point<2>> inner_border)
    {
        CDT cdt;
        create_grid_using_cgal<cguc_type::Default>(
                cdt, outer_border, inner_border);

        // std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        // std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
        // 
        // size_t mesh_faces_counter = 0;

        // // for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
        // for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
        // {
        //     // if(fit->is_in_domain()) 
        //     //     ++mesh_faces_counter;
        // };
        // 
        // std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;

        convert_to_dealii_format(cdt, triangulation);
    };

    void set_grid(
            dealii::Triangulation< 2 > &triangulation,
            vec<prmt::Point<2>> outer_border,
            vec<prmt::Point<2>> inner_border,
            vec<st> type_outer_border)
    {
        CDT cdt;
        vec<CDT::Point> border_on_cdt;
        create_grid_using_cgal_1(cdt, border_on_cdt,
                outer_border, inner_border);

        // std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
        // std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
        // 
        // size_t mesh_faces_counter = 0;

        // // for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
        // for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
        // {
        //     // if (fit->is_in_domain()) 
        //     //     ++mesh_faces_counter;
        // };
        // 
        // std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;

        convert_to_dealii_format(cdt, border_on_cdt,
                triangulation, outer_border, type_outer_border);
    };

    void set_grid(
            dealii::Triangulation< 2 > &triangulation,
            vec<prmt::Point<2>> outer_border,
            vec<vec<prmt::Point<2>>> inner_border,
            vec<st> type_outer_border)
    {
        CDT cdt;
        vec<CDT::Point> border_on_cdt;
        create_grid_using_cgal_1(cdt, border_on_cdt,
                outer_border, inner_border);

        convert_to_dealii_format(cdt, border_on_cdt,
                triangulation, outer_border, type_outer_border);
    };

    void make_grid(
            dealii::Triangulation< 2 > &triangulation,
            vec<prmt::Point<2>> outer_border,
            vec<st> type_outer_border)
    {
        CDT cdt;
        vec<CDT::Point> border_on_cdt;
        printf("2.1\n");
        create_grid_using_cgal_2(cdt, border_on_cdt, outer_border);
        printf("2.2\n");

        convert_to_dealii_format_without_domains(cdt, border_on_cdt,
                triangulation, outer_border, type_outer_border);
        printf("2.3\n");
    };

};

#endif
