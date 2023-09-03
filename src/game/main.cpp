#include "gameutils/edge_soup.h"

const ivec2 screen_size = ivec2(480, 270);
const std::string_view window_name = "2D contours";

Interface::Window window(std::string(window_name), screen_size * 2, Interface::windowed, adjust_(Interface::WindowSettings{}, .min_size = screen_size));
static Graphics::DummyVertexArray dummy_vao = nullptr;

Audio::Context audio_context = nullptr;
Audio::SourceManager audio_controller;

const Graphics::ShaderConfig shader_config = Graphics::ShaderConfig::Core();
Interface::ImGuiController gui_controller(Poly::derived<Interface::ImGuiController::GraphicsBackend_Modern>, adjust_(Interface::ImGuiController::Config{}, .shader_header = shader_config.common_header, .store_state_in_file = {}));

namespace Fonts
{
    namespace Files
    {
        Graphics::FontFile main(Program::ExeDir() + "assets/Monocat_6x12.ttf", 12);
    }
    Graphics::Font main;
}

GameUtils::AdaptiveViewport adaptive_viewport(shader_config, screen_size);
Render r = adjust_(Render(0x2000, shader_config), .SetMatrix(adaptive_viewport.GetDetails().MatrixCentered()));

Input::Mouse mouse;

Random::DefaultGenerator random_generator = Random::MakeGeneratorFromRandomDevice();
Random::DefaultInterfaces<Random::DefaultGenerator> ra(random_generator);

template class EdgeSoup<int>;
template class EdgeSoup<float>;

struct ContourDemo
{
    using scalar = float;
    using vector = vec2<scalar>;
    using matrix = mat2<scalar>;

    std::vector<vector> unfinished_contour;

    using Shape = EdgeSoup<scalar>;
    Shape shape;
    vector self_pos;
    matrix self_matrix;

    Shape other_shape;
    vector other_pos;
    matrix other_matrix;

    void DrawLine(fvec2 a, fvec2 b, fvec3 color, float alpha = 1) const
    {
        r.fquad(a, fvec2((b - a).len(), 1)).center(fvec2(0,0.5f)).color(color).alpha(alpha).rotate((b - a).angle());
    }
    void DrawLineWithNormal(fvec2 a, fvec2 b, fvec3 color, float alpha = 1) const
    {
        DrawLine(a, b, color, alpha);
        fvec2 center = (a + b) / 2;
        DrawLine(center, center + (b - a).rot90(-1).norm() * 8, color, alpha);
    }
    void DrawRect(frect2 rect, fvec2 pos, fmat2 matrix, fvec3 color, float alpha = 1) const
    {
        vector prev = matrix * rect.corner(3) + pos;
        for (vector point : rect.to_contour())
        {
            point = matrix * point + pos;
            DrawLine(prev, point, color, alpha);
            prev = point;
        }
    }
    void DrawShape(bool draw_aabb, const Shape &target_shape, Shape::vector target_pos, Shape::matrix target_rot, fvec3 color, float alpha = 1) const
    {
        auto bounds = target_shape.Bounds().to_contour();
        for (auto &vertex : bounds)
            vertex = target_rot * vertex + target_pos;

        if (draw_aabb)
        {
            r.ftriangle(bounds[0], bounds[1], bounds[3]).color(fvec3(1,1,1)).alpha(0.1f);
            r.ftriangle(bounds[1], bounds[2], bounds[3]).color(fvec3(1,1,1)).alpha(0.1f);
        }

        target_shape.EdgeTree().CollideCustom([&](auto &&){return true;}, [&](Shape::AabbTree::NodeIndex e)
        {
            Shape::Edge edge = target_shape.GetEdge(Shape::EdgeIndex(e));
            edge.a = target_rot * edge.a;
            edge.b = target_rot * edge.b;
            edge.a += target_pos;
            edge.b += target_pos;
            DrawLineWithNormal(edge.a, edge.b, color, alpha);
            // Edge IDs.
            if (false)
            {
                auto maybe_round = []<typename T>(vec2<T> value)
                {
                    if constexpr (Math::floating_point_scalar<T>)
                        return iround(value);
                    else
                        return value;
                };
                auto edge_center = maybe_round((edge.a + edge.b) / 2);
                for (int i = 0; i < 4; i++)
                    r.itext(edge_center + ivec2::dir4(i), Graphics::Text(Fonts::main, FMT("{}", e))).align(ivec2(0)).color(fvec3(0,0,0));
                r.itext(edge_center, Graphics::Text(Fonts::main, FMT("{}", e))).align(ivec2(0)).color(fvec3(1,1,1));
            }
            return false;
        });
    }

    void AddUnfinishedContourToShape(Shape *target_shape = nullptr)
    {
        if (!target_shape)
            target_shape = &shape;

        if (!unfinished_contour.empty())
        {
            std::optional<vector> prev;
            for (vector point : unfinished_contour)
            {
                if (prev)
                    target_shape->AddEdge({.a = *prev, .b = point});
                prev = point;
            }
            target_shape->AddEdge({.a = prev.value(), .b = unfinished_contour.front()});

            std::cout << Refl::ToString(unfinished_contour) << '\n';

            unfinished_contour.clear();
        }
    }

    void Init()
    {
        self_pos = vector(120, 80);
        self_matrix = matrix::rotate(to_rad(-60));

        other_pos = vector(-9, 77);
        // other_matrix = imat2(0,-1,1,0);
        other_matrix = matrix::rotate(to_rad(10));

        // C
        // Refl::FromString(unfinished_contour, "[(-70,-65),(21,-81),(35,6),(-58,21),(-44,-4),(10,-12),(1,-59),(-56,-49)]");
        // flipped C
        // Refl::FromString(unfinished_contour, "[(56,-49),(-1,-59),(-10,-12),(44,-4),(58,21),(-35,6),(-21,-81),(70,-65)]");
        // O
        // Refl::FromString(unfinished_contour, "[(-60,-50),(-11,-68),(41,-48),(58,6),(37,51),(-15,67),(-63,49),(-85,0)]");
        // overhangs
        // Refl::FromString(unfinished_contour, "[(-26,-121),(12,-46),(-9,-40),(10,-6),(-11,-6),(9,28),(-8,23),(7,56),(-6,56),(-17,31),(-34,119),(-100,128),(-102,-128)]");
        // hook
        // Refl::FromString(unfinished_contour, "[(-67,31),(-67,14),(-34,25),(-8,10),(-8,-37),(-17,-49),(1,-38),(10,20),(-27,48)]");
        // rod
        // Refl::FromString(unfinished_contour, "[(-3,22),(-24,29),(-33,10),(-14,-124),(-5,-118)]");

        // flipped C
        Refl::FromString(unfinished_contour, "[(56,-49),(-1,-59),(-10,-12),(44,-4),(58,21),(-35,6),(-21,-81),(70,-65)]");
        AddUnfinishedContourToShape();

        Refl::FromString(unfinished_contour, "[(-3,22),(-24,29),(-5,-118)]"); // rod
        AddUnfinishedContourToShape(&other_shape);
    }

    void Tick()
    {
        if (mouse.left.pressed())
            unfinished_contour.push_back(mouse.pos());

        if (mouse.right.pressed())
            AddUnfinishedContourToShape();

        other_pos = mouse.pos();
    }

    void Render() const
    {
        // Test.
        r.itext(-screen_size/2, Graphics::Text(Fonts::main, FMT("cursor = {}", mouse.pos()))).align(ivec2(-1)).color(fvec3(1));

        { // The shape.
            bool collides = false;
            vector self_offset;
            matrix self_matrix_offset;
            vector other_offset;
            matrix other_matrix_offset;
            if (unfinished_contour.empty())
            {
                // // Point collision.
                // auto collider = shape.MakePointCollider(vector(0));
                // collides = collider.CollidePoint(mouse.pos());
                // r.fquad(mouse.pos() + 0.5f, fvec2(32, 1)).center(fvec2(0.5f)).rotate(collider.DebugRayDirection() * f_pi / 2).color(fvec3(0,0,1));

                // // Shape collision.
                // collides = shape.CollideEdgeSoupSimple(other_shape, other_pos, other_matrix);

                // Shape collision with offset.
                int num_steps = 8;
                vector self_vel(-31, 12);
                vector other_vel(40, 20);
                float self_ang_vel = -f_pi/4;
                float other_ang_vel = f_pi;
                self_matrix_offset = matrix::rotate(self_ang_vel);
                other_matrix_offset = matrix::rotate(other_ang_vel);

                self_offset = self_vel;
                other_offset = other_vel;

                Storage::MonotonicPool persistent_pool;
                Storage::MonotonicPool temp_pool;

                auto collider = shape.MakeEdgeSoupCollider(other_shape, {
                    .persistent_pool = &persistent_pool,
                    .temp_pool = &temp_pool,
                    .self_pos = self_pos,
                    .other_pos = other_pos,
                    .self_vel = self_vel,
                    .other_vel = other_vel,
                    .self_rot = self_matrix,
                    .other_rot = other_matrix,
                    .self_angular_vel_abs_upper_bound = abs(self_ang_vel),
                    .other_angular_vel_abs_upper_bound = abs(other_ang_vel),
                });

                // collides = collider.Collide(vector{}, matrix{}, other_pos, other_matrix);

                if (!collider.Collide(self_pos + self_offset, self_matrix_offset * self_matrix, other_pos + other_offset, other_matrix_offset * other_matrix))
                {
                    collides = false;
                }
                else
                {
                    float t = 0.5f;
                    self_offset = self_vel * t;
                    self_matrix_offset = matrix::rotate(self_ang_vel * t);
                    other_offset = other_vel * t;
                    other_matrix_offset = matrix::rotate(other_ang_vel * t);
                    for (int i = 0; i < num_steps; i++)
                    {
                        bool hit = collider.Collide(self_pos + self_offset, self_matrix_offset * self_matrix, other_pos + other_offset, other_matrix_offset * other_matrix);
                        if (hit)
                            collides = true;

                        t += (hit ? -1 : 1) * (1.f / (1 << (i + 2)));
                        self_offset = self_vel * t;
                        other_offset = other_vel * t;
                        self_matrix_offset = matrix::rotate(self_ang_vel * t);
                        other_matrix_offset = matrix::rotate(other_ang_vel * t);
                    }
                }
            }

            DrawShape(false, shape, self_pos + self_offset, self_matrix_offset * self_matrix, fvec3(1,1,1), 0.5f);
            DrawShape(true, shape, self_pos, self_matrix, collides ? fvec3(1,0,1) : fvec3(1,1,1));

            DrawShape(false, other_shape, other_pos + other_offset, other_matrix_offset * other_matrix, fvec3(0,0.5f,1), 0.5f);
            DrawShape(true, other_shape, other_pos, other_matrix, fvec3(0,0.5f,1));
        }

        { // Unfinished contour.
            std::optional<vector> prev;
            for (vector point : unfinished_contour)
            {
                if (prev)
                    DrawLineWithNormal(*prev, point, fvec3(0,1,0));
                prev = point;
            }
            if (prev)
            {
                DrawLineWithNormal(*prev, mouse.pos(), mouse.pos() - *prev == any(0) ? fvec3(1,0.5f,0) : fvec3(1,1,0));
                DrawLineWithNormal(*prev, unfinished_contour.front(), fvec3(1,1,1), 0.3f);
            }
        }
    }
};

struct Application : Program::DefaultBasicState
{
    GameUtils::FpsCounter fps_counter;

    ContourDemo demo;

    void Resize()
    {
        adaptive_viewport.Update();
        mouse.SetMatrix(adaptive_viewport.GetDetails().MouseMatrixCentered());
    }

    Metronome metronome = Metronome(60);

    Metronome *GetTickMetronome() override
    {
        return &metronome;
    }

    int GetFpsCap() override
    {
        return 60 * NeedFpsCap();
    }

    void EndFrame() override
    {
        fps_counter.Update();
        window.SetTitle(STR((window_name), " TPS:", (fps_counter.Tps()), " FPS:", (fps_counter.Fps())));
    }

    void Tick() override
    {
        // window.ProcessEvents();
        window.ProcessEvents({gui_controller.EventHook()});

        if (window.ExitRequested())
            Program::Exit();
        if (window.Resized())
        {
            Resize();
            Graphics::Viewport(window.Size());
        }

        gui_controller.PreTick();
        demo.Tick();
        audio_controller.Tick();

        Audio::CheckErrors();
    }

    void Render() override
    {
        gui_controller.PreRender();
        adaptive_viewport.BeginFrame();
        Graphics::SetClearColor(fvec3(0));
        Graphics::Clear();
        r.BindShader();
        demo.Render();
        r.Finish();
        adaptive_viewport.FinishFrame();
        gui_controller.PostRender();
        Graphics::CheckErrors();

        window.SwapBuffers();
    }


    void Init()
    {
        ImGui::StyleColorsDark();

        { // Load images.
            Graphics::GlobalData::Load(Program::ExeDir() + "assets/images/");
            r.SetAtlas("");

            // Load the font atlas.
            struct FontLoader : Graphics::GlobalData::Generator
            {
                ivec2 Size() const override
                {
                    return ivec2(256);
                }

                void Generate(Graphics::Image &image, irect2 rect) override
                {
                    Unicode::CharSet glyph_ranges;
                    glyph_ranges.Add(Unicode::Ranges::Basic_Latin);

                    Graphics::MakeFontAtlas(image, rect, {
                        {Fonts::main, Fonts::Files::main, glyph_ranges, Graphics::FontFile::monochrome_with_hinting},
                    });
                }
            };
            (void)Graphics::GlobalData::Image<"font_atlas", FontLoader>();
        }

        // Load various small fonts
        auto monochrome_font_flags = ImGuiFreeTypeBuilderFlags_Monochrome | ImGuiFreeTypeBuilderFlags_MonoHinting;

        gui_controller.LoadFont(Program::ExeDir() + "assets/Monocat_6x12.ttf", 12.0f, adjust(ImFontConfig{}, .FontBuilderFlags = monochrome_font_flags));
        gui_controller.LoadDefaultFont();

        Graphics::Blending::Enable();
        Graphics::Blending::FuncNormalPre();

        demo.Init();
    }
};

IMP_MAIN(,)
{
    Application app;
    app.Init();
    app.Resize();
    app.RunMainLoop();
    return 0;
}
