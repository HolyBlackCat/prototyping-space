#pragma once

#include <cstddef>
#include <utility>

#include <cglfl/cglfl.hpp>

#include "graphics/image.h"
#include "macros/finally.h"
#include "utils/mat.h"
#include "utils/sparse_set.h"

namespace Graphics
{
    enum InterpolationMode
    {
        nearest,
        linear,
        min_nearest_mag_linear,
        min_linear_mag_nearest,
    };

    enum WrapMode
    {
        clamp  = GL_CLAMP_TO_EDGE,
        mirror = GL_MIRRORED_REPEAT,
        repeat = GL_REPEAT,
        #ifdef GL_CLAMP_TO_BORDER
        fill   = GL_CLAMP_TO_BORDER,
        #endif
    };

    class TexObject
    {
        struct Data
        {
            GLuint handle = 0;
        };

        Data data;

      public:
        TexObject() {}

        TexObject(decltype(nullptr))
        {
            glGenTextures(1, &data.handle);
            if (!data.handle)
                throw std::runtime_error("Unable to create a texture.");
            // Not needed because there is no code below this point:
            // FINALLY_ON_THROW{glDeleteTextures(1, &data.handle);};
        }

        TexObject(TexObject &&other) noexcept : data(std::exchange(other.data, {})) {}
        TexObject &operator=(TexObject other) noexcept
        {
            std::swap(data, other.data);
            return *this;
        }

        ~TexObject()
        {
            // Deleting a texture unbinds it.
            if (data.handle)
                glDeleteTextures(1, &data.handle); // Deleting 0 is a no-op, but GL could be unloaded at this point.
        }

        explicit operator bool() const
        {
            return bool(data.handle);
        }

        GLuint Handle() const
        {
            return data.handle;
        }
    };

    class TexUnit
    {
        using index_alloc_t = SparseSet<int>;

        static index_alloc_t &Allocator() // Wrapped into a function to prevent the static init order fiasco.
        {
            static index_alloc_t ret(64);
            return ret;
        }

        struct Data
        {
            int index = -1;
            GLuint handle = 0;
        };

        Data data;

        inline static int active_index = 0;

      public:
        TexUnit()
        {
            // We need to create the allocator even in the null constructor to dodge the destruction order fiasco.
            (void)Allocator();
        }

        TexUnit(decltype(nullptr))
        {
            if (Allocator().IsFull())
                throw std::runtime_error("No free texture units.");
            data.index = Allocator().InsertAny();
        }
        explicit TexUnit(GLuint handle) : TexUnit(nullptr)
        {
            AttachHandle(handle);
        }
        explicit TexUnit(const TexObject &texture) : TexUnit(nullptr)
        {
            Attach(texture);
        }

        TexUnit(TexUnit &&other) noexcept : data(std::exchange(other.data, {})) {}
        TexUnit &operator=(TexUnit other) noexcept
        {
            std::swap(data, other.data);
            return *this;
        }

        ~TexUnit()
        {
            if (*this)
                Allocator().EraseUnordered(data.index);
        }

        explicit operator bool() const
        {
            return data.index != -1;
        }

        int Index() const
        {
            return data.index;
        }

        static void ActivateIndex(int index)
        {
            if (active_index == index)
                return;
            glActiveTexture(GL_TEXTURE0 + index);
            active_index = index;
        }
        void Activate()
        {
            ASSERT(*this, "Attempt to use a null texture unit.");
            if (!*this)
                return;
            ActivateIndex(data.index);
        }
        [[nodiscard]] bool Active()
        {
            return bool(*this) && active_index == data.index;
        }

        TexUnit &&AttachHandle(GLuint handle)
        {
            ASSERT(*this, "Attempt to use a null texture unit.");
            if (!*this)
                return std::move(*this);

            Activate();
            data.handle = handle;
            glBindTexture(GL_TEXTURE_2D, handle);
            return std::move(*this);
        }
        TexUnit &&Attach(const TexObject &texture)
        {
            AttachHandle(texture.Handle());
            return std::move(*this);
        }
        TexUnit &&Detach()
        {
            AttachHandle(0);
            return std::move(*this);
        }

        GLuint AttachedHandle() const
        {
            return data.handle;
        }
        bool HasAttachedHandle() const
        {
            return bool(data.handle);
        }

        TexUnit &&Interpolation(InterpolationMode mode)
        {
            ASSERT(HasAttachedHandle(), "Attempt to use a texture unit without an attached texture.");
            if (!HasAttachedHandle())
                return std::move(*this);

            Activate();

            GLenum min_mode = (mode == nearest || mode == min_nearest_mag_linear ? GL_NEAREST : GL_LINEAR);
            GLenum mag_mode = (mode == nearest || mode == min_linear_mag_nearest ? GL_NEAREST : GL_LINEAR);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_mode);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_mode);

            return std::move(*this);
        }

        TexUnit &&WrapX(WrapMode mode)
        {
            ASSERT(HasAttachedHandle(), "Attempt to use a texture unit without an attached texture.");
            if (!HasAttachedHandle())
                return std::move(*this);

            Activate();
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GLuint(mode));
            return std::move(*this);
        }
        TexUnit &&WrapY(WrapMode mode)
        {
            ASSERT(HasAttachedHandle(), "Attempt to use a texture unit without an attached texture.");
            if (!HasAttachedHandle())
                return std::move(*this);

            Activate();
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GLuint(mode));
            return std::move(*this);
        }
        TexUnit &&Wrap(WrapMode mode)
        {
            WrapX(mode);
            WrapY(mode);
            return std::move(*this);
        }

        TexUnit &&SetData(ivec2 size, const uint8_t *pixels = 0)
        {
            GLenum internal_format =
            #ifdef GL_RGBA8
                GL_RGBA8;
            #else
                GL_RGBA;
            #endif

            SetData(internal_format, GL_RGBA, GL_UNSIGNED_BYTE, size, pixels);
            return std::move(*this);
        }
        TexUnit &&SetData(GLenum internal_format, GLenum format, GLenum type, ivec2 size, const uint8_t *pixels = 0)
        {
            ASSERT(HasAttachedHandle(), "Attempt to use a texture unit without an attached texture.");
            if (!HasAttachedHandle())
                return std::move(*this);

            Activate();
            glTexImage2D(GL_TEXTURE_2D, 0, internal_format, size.x, size.y, 0, format, type, pixels);
            return std::move(*this);
        }
        TexUnit &&SetData(const Image &image)
        {
            ASSERT(image, "Attempt to use a null image.");
            SetData(image.Size(), image.Data());
            return std::move(*this);
        }

        TexUnit &&SetDataPart(ivec2 pos, ivec2 size, const uint8_t *pixels)
        {
            SetDataPart(GL_RGBA, GL_UNSIGNED_BYTE, pos, size, pixels);
            return std::move(*this);
        }
        TexUnit &&SetDataPart(GLenum format, GLenum type, ivec2 pos, ivec2 size, const uint8_t *pixels)
        {
            ASSERT(HasAttachedHandle(), "Attempt to use a texture unit without an attached texture.");
            if (!HasAttachedHandle())
                return std::move(*this);

            Activate();
            glTexSubImage2D(GL_TEXTURE_2D, 0, pos.x, pos.y, size.x, size.y, format, type, pixels);
            return std::move(*this);
        }
    };

    class Texture
    {
        TexObject object;
        TexUnit unit;
        ivec2 size = ivec2(0);

      public:
        Texture() {}

        Texture(decltype(nullptr)) : object(nullptr), unit(nullptr)
        {
            unit.Attach(object);
        }

        Texture &&Interpolation(InterpolationMode mode)
        {
            unit.Interpolation(mode);
            return std::move(*this);
        }

        Texture &&WrapX(WrapMode mode)
        {
            unit.WrapX(mode);
            return std::move(*this);
        }
        Texture &&WrapY(WrapMode mode)
        {
            unit.WrapY(mode);
            return std::move(*this);
        }
        Texture &&Wrap(WrapMode mode)
        {
            unit.Wrap(mode);
            return std::move(*this);
        }

        Texture &&SetData(ivec2 new_size, const uint8_t *pixels = 0)
        {
            unit.SetData(new_size, pixels);
            size = new_size;
            return std::move(*this);
        }
        Texture &&SetData(GLenum internal_format, GLenum format, GLenum type, ivec2 new_size, const uint8_t *pixels = 0)
        {
            unit.SetData(internal_format, format, type, new_size, pixels);
            size = new_size;
            return std::move(*this);
        }
        Texture &&SetData(const Image &image)
        {
            unit.SetData(image);
            size = image.Size();
            return std::move(*this);
        }

        Texture &&SetDataPart(ivec2 part_pos, ivec2 part_size, const uint8_t *pixels)
        {
            unit.SetDataPart(part_pos, part_size, pixels);
            return std::move(*this);
        }
        Texture &&SetDataPart(GLenum format, GLenum type, ivec2 part_pos, ivec2 part_size, const uint8_t *pixels)
        {
            unit.SetDataPart(format, type, part_pos, part_size, pixels);
            return std::move(*this);
        }

        ivec2 Size() const
        {
            return size;
        }
        int Handle() const
        {
            return object.Handle();
        }
        int Index() const
        {
            return unit.Index();
        }

              TexObject &Object()       {return object;}
        const TexObject &Object() const {return object;}

              TexUnit &Unit()       {return unit;}
        const TexUnit &Unit() const {return unit;}

        operator       TexObject &()       {return object;}
        operator const TexObject &() const {return object;}

        operator       TexUnit &()       {return unit;}
        operator const TexUnit &() const {return unit;}
    };
}
