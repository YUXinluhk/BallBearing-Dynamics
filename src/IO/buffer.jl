# =====================================================================
# IO/buffer.jl — Ring buffer for async HDF5 field output persistence
#
# Prevents OOM by flushing accumulated field output to disk periodically.
# =====================================================================

using HDF5

"""
    OutputBuffer — accumulates field output frames, flushes to HDF5.
"""
mutable struct OutputBuffer
    data::Matrix{Float64}    # (n_fields, buffer_capacity)
    write_idx::Int           # next write position (1-indexed)
    total_written::Int       # total frames flushed to disk
    filepath::String
    dataset_name::String
    n_fields::Int
    capacity::Int
end

"""
    OutputBuffer(filepath, n_fields; capacity=1000)
"""
function OutputBuffer(filepath::String, n_fields::Int;
                      capacity::Int=1000, dataset_name::String="field_output")
    OutputBuffer(
        zeros(n_fields, capacity),
        1, 0,
        filepath, dataset_name,
        n_fields, capacity,
    )
end

"""
    push_frame!(buf, frame)

Add one field output frame. Auto-flushes when buffer is full.
"""
function push_frame!(buf::OutputBuffer, frame::AbstractVector{Float64})
    @assert length(frame) == buf.n_fields "Frame length mismatch"

    buf.data[:, buf.write_idx] .= frame
    buf.write_idx += 1

    if buf.write_idx > buf.capacity
        flush_buffer!(buf)
    end

    return nothing
end

"""
    flush_buffer!(buf)

Write buffered data to HDF5 (append or create).
"""
function flush_buffer!(buf::OutputBuffer)
    n_frames = buf.write_idx - 1
    n_frames <= 0 && return nothing

    chunk = buf.data[:, 1:n_frames]

    h5open(buf.filepath, isfile(buf.filepath) ? "r+" : "w") do fid
        ds_name = buf.dataset_name

        if haskey(fid, ds_name)
            # Extend existing dataset
            ds = fid[ds_name]
            existing_size = size(ds, 2)
            new_size = existing_size + n_frames
            HDF5.set_extent_dims(ds, (buf.n_fields, new_size))
            ds[:, existing_size+1:new_size] = chunk
        else
            # Create chunked, resizable dataset
            ds = create_dataset(fid, ds_name, Float64,
                                (buf.n_fields, n_frames),
                                max_dims=(buf.n_fields, -1),
                                chunk=(buf.n_fields, min(n_frames, 100)))
            ds[:, :] = chunk
        end
    end

    buf.total_written += n_frames
    buf.write_idx = 1  # reset pointer

    return nothing
end

"""
    finalize_buffer!(buf)

Flush remaining data and close.
"""
function finalize_buffer!(buf::OutputBuffer)
    flush_buffer!(buf)
    return buf.total_written
end
