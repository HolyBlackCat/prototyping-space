#include "archive.h"

#include <type_traits>

#include <zlib.h>

#include "utils/robust_math.h"

namespace Archive
{
    namespace Raw
    {
        std::size_t MaxCompressedSize(const uint8_t *src_begin, const uint8_t *src_end)
        {
            return compressBound(src_end - src_begin);
        }

        uint8_t *Compress(const uint8_t *src_begin, const uint8_t *src_end, uint8_t *dst_begin, uint8_t *dst_end)
        {
            uLong dst_size = dst_end - dst_begin; // compress() changes this value.
            int status = compress(dst_begin, &dst_size, src_begin, src_end - src_begin);
            if (status != Z_OK)
                throw std::runtime_error("Compression failure.");
            return dst_begin + dst_size;
        }

        void Uncompress(const uint8_t *src_begin, const uint8_t *src_end, uint8_t *dst_begin, uint8_t *dst_end)
        {
            uLong dst_size = dst_end - dst_begin; // uncompress() changes this value.
            int status = uncompress(dst_begin, &dst_size, src_begin, src_end - src_begin);
            if (status != Z_OK || dst_size != uLong(dst_end - dst_begin))
                throw std::runtime_error("Uncompression failure.");
        }
    }


    using size_type = uint64_t;
    static_assert(std::is_unsigned_v<size_type> && sizeof(size_type) >= sizeof(std::size_t), "`size_type` must be an unsigned type not smaller than `std::size_t`.");

    [[nodiscard]] std::size_t MaxCompressedSize(const uint8_t *src_begin, const uint8_t *src_end)
    {
        return sizeof(size_type) + Raw::MaxCompressedSize(src_begin, src_end);
    }

    [[nodiscard]] uint8_t *Compress(const uint8_t *src_begin, const uint8_t *src_end, uint8_t *dst_begin, uint8_t *dst_end)
    {
        if (dst_end - dst_begin < std::ptrdiff_t(sizeof(size_type)))
            throw std::runtime_error("Compression failure.");

        std::size_t size = src_end - src_begin;

        if (size > std::numeric_limits<size_type>::max())
            throw std::runtime_error("Compression failure.");

        for (std::size_t i = 0; i < sizeof(size_type); i++)
            dst_begin[i] = (size >> (i * 8)) & 0xff;

        return Raw::Compress(src_begin, src_end, dst_begin + sizeof(size_type), dst_end);
    }

    [[nodiscard]] std::size_t UncompressedSize(const uint8_t *src_begin, const uint8_t *src_end)
    {
        if (src_end - src_begin < std::ptrdiff_t(sizeof(size_type)))
            throw std::runtime_error("Uncompression failure.");

        size_type size = 0;
        for (std::size_t i = 0; i < sizeof(size_type); i++)
            size |= (size_type(src_begin[i]) << (i * 8));

        std::size_t ret;
        if (Robust::conversion_fails(size, ret))
            throw std::runtime_error("Unable to uncompress: The object is too large.");

        return ret;
    }

    void Uncompress(const uint8_t *src_begin, const uint8_t *src_end, uint8_t *dst_begin)
    {
        std::size_t size = UncompressedSize(src_begin, src_end);
        Raw::Uncompress(src_begin + sizeof(size_type), src_end, dst_begin, dst_begin + size);
    }
}
